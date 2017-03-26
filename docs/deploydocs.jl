using ..Documenter:
        Utilities,
        Writers,
        Deps
function Documenter.deploydocs(;
        root   = Utilities.currentdir(),
        target = "site",
        dirname = "",

        repo   = error("no 'repo' keyword provided."),
        branch = "gh-pages",
        latest = "master",

        osname = "linux",
        julia  = "nightly",

        deps   = Deps.pip("pygments", "mkdocs"),
        make   = () -> run(`mkdocs build`),

        figs_dir = "assets/figures",
    )
    # Get environment variables.
    documenter_key      = get(ENV, "DOCUMENTER_KEY",       "")
    travis_branch       = get(ENV, "TRAVIS_BRANCH",        "")
    travis_pull_request = get(ENV, "TRAVIS_PULL_REQUEST",  "")
    travis_repo_slug    = get(ENV, "TRAVIS_REPO_SLUG",     "")
    travis_tag          = get(ENV, "TRAVIS_TAG",           "")
    travis_osname       = get(ENV, "TRAVIS_OS_NAME",       "")
    travis_julia        = get(ENV, "TRAVIS_JULIA_VERSION", "")

    # Should figures be re-drawn (default to true)'
    draw_fig            = get(ENV, "DRAW_FIG",             "")
    draw_fig = !(draw_fig == "false")

    # Other variables.
    sha = cd(root) do
        # We'll make sure we run the git commands in the source directory (root), in case
        # the working directory has been changed (e.g. if the makedocs' build argument is
        # outside root).
        try
            readchomp(`git rev-parse --short HEAD`)
        catch
            # git rev-parse will throw an error and return code 128 if it is not being
            # run in a git repository, which will make run/readchomp throw an exception.
            # We'll assume that if readchomp fails it is due to this and set the sha
            # variable accordingly.
            "(not-git-repo)"
        end
    end

    # Sanity checks
    if !isa(julia, AbstractString)
        error("julia must be a string, got $julia ($(typeof(julia)))")
    end
    if !isempty(travis_repo_slug) && !contains(repo, travis_repo_slug)
        warn("repo $repo does not match $travis_repo_slug")
    end

    # When should a deploy be attempted?
    should_deploy =
        contains(repo, travis_repo_slug) &&
        travis_pull_request == "false"   &&
        travis_osname == osname &&
        travis_julia  == julia  &&
        (
            travis_branch == latest ||
            travis_tag    != ""
        )

    # check DOCUMENTER_KEY only if the branch, Julia version etc. check out
    if should_deploy && isempty(documenter_key)
        warn("""
            DOCUMENTER_KEY environment variable missing, unable to deploy.
              Note that in Documenter v0.9.0 old deprecated authentication methods were removed.
              DOCUMENTER_KEY is now the only option. See the documentation for more information.""")
        should_deploy = false
    end

    if get(ENV, "DOCUMENTER_DEBUG", "") == "true"
        Utilities.debug("TRAVIS_REPO_SLUG       = \"$travis_repo_slug\"")
        Utilities.debug("  should match \"$repo\" (kwarg: repo)")
        Utilities.debug("TRAVIS_PULL_REQUEST    = \"$travis_pull_request\"")
        Utilities.debug("  deploying if equal to \"false\"")
        Utilities.debug("TRAVIS_OS_NAME         = \"$travis_osname\"")
        Utilities.debug("  deploying if equal to \"$osname\" (kwarg: osname)")
        Utilities.debug("TRAVIS_JULIA_VERSION   = \"$travis_julia\"")
        Utilities.debug("  deploying if equal to \"$julia\" (kwarg: julia)")
        Utilities.debug("TRAVIS_BRANCH          = \"$travis_branch\"")
        Utilities.debug("TRAVIS_TAG             = \"$travis_tag\"")
        Utilities.debug("  deploying if branch equal to \"$latest\" (kwarg: latest) or tag is set")
        Utilities.debug("git commit SHA         = $sha")
        Utilities.debug("DOCUMENTER_KEY exists  = $(!isempty(documenter_key))")
        Utilities.debug("should_deploy          = $should_deploy")
    end

    if should_deploy
        # Add local bin path if needed.
        Deps.updatepath!()
        # Install dependencies when applicable.
        if deps !== nothing
            Utilities.log("installing dependencies.")
            deps()
        end
        # Change to the root directory and try to deploy the docs.
        cd(root) do
            Utilities.log("setting up target directory.")
            isdir(target) || mkpath(target)
            # Run extra build steps defined in `make` if required.
            if make !== nothing
                Utilities.log("running extra build steps.")
                make()
            end
            Utilities.log("pushing new documentation to remote: $repo:$branch.")
            mktempdir() do temp
                dirname = isempty(dirname) ? temp : joinpath(temp, dirname)
                isdir(dirname) || mkpath(dirname)
                # Versioned docs directories.
                latest_dir = joinpath(dirname, "latest")
                stable_dir = joinpath(dirname, "stable")
                tagged_dir = joinpath(dirname, travis_tag)

                keyfile = abspath(joinpath(root, ".documenter"))
                target_dir = abspath(target)

                # The upstream URL to which we push new content and the ssh decryption commands.
                write(keyfile, Compat.String(base64decode(documenter_key)))
                chmod(keyfile, 0o600)
                upstream = "git@$(replace(repo, "github.com/", "github.com:"))"

                # Use a custom SSH config file to avoid overwriting the default user config.
                Documenter.withfile(joinpath(homedir(), ".ssh", "config"),
                    """
                    Host github.com
                        StrictHostKeyChecking no
                        HostName github.com
                        IdentityFile $keyfile
                    """
                ) do
                    cd(temp) do
                        # Setup git.
                        run(`git init`)
                        run(`git config user.name "autodocs"`)
                        run(`git config user.email "autodocs"`)

                        # Fetch from remote and checkout the branch.
                        success(`git remote add upstream $upstream`) ||
                            error("could not add new remote repo.")

                        success(`git fetch upstream`) ||
                            error("could not fetch from remote.")

                        success(`git checkout -b $branch upstream/$branch`) ||
                            error("could not checkout remote branch.")

                        # Copy docs to `latest`, or `stable`, `<release>`, and `<version>` directories.
                        if isempty(travis_tag)
                            if !draw_fig
                                # If figures not drawn, then keep old ones
                                # (only for latest: always draw figures for stable/tagged)
                                cp(latest_dir*figs_dir, target_dir*figs_dir; remove_destination = true)
                            end
                            cp(target_dir, latest_dir; remove_destination = true)
                        else
                            cp(target_dir, stable_dir; remove_destination = true)
                            cp(target_dir, tagged_dir; remove_destination = true)
                            # Build a `release-*.*` folder as well when the travis tag is
                            # valid, which it *should* always be anyway.
                            if ismatch(Base.VERSION_REGEX, travis_tag)
                                local version = VersionNumber(travis_tag)
                                local release = "release-$(version.major).$(version.minor)"
                                cp(target_dir, joinpath(dirname, release); remove_destination = true)
                            end
                        end

                        # Create the versions.js file containing a list of all docs
                        # versions. This must always happen after the folder copying.
                        Writers.HTMLWriter.generate_version_file(dirname)

                        # Add, commit, and push the docs to the remote.
                        run(`git add -A .`)
                        try run(`git commit -m "build based on $sha"`) end

                        success(`git push -q upstream HEAD:$branch`) ||
                            error("could not push to remote repo.")

                        # Remove the unencrypted private key.
                        isfile(keyfile) && rm(keyfile)
                    end
                end
            end
        end
    else
        Utilities.log("""
            skipping docs deployment.
              You can set DOCUMENTER_DEBUG to "true" in Travis to see more information.""")
    end
end
