#benchmarking in development
# see https://juliaci.github.io/PkgBenchmark.jl/stable/run_benchmarks.html
# https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/doc/manual.md
# seehttps://github.com/JuliaArrays/StaticArrays.jl for good examples

# To upload to GitHub as markdown: #TODO 
#https://juliaci.github.io/PkgBenchmark.jl/stable/export_markdown.html
using BenchmarkTools

#to run:
# using PkgBenchmark
# benchmarkpkg("PhyloNetworks")

# Define a parent BenchmarkGroup to contain our SUITE
SUITE = BenchmarkGroup()

# EXAMPLES
# Add some child groups to our benchmark SUITE. The most relevant BenchmarkGroup constructor
# for this case is BenchmarkGroup(tags::Vector). These tags are useful for
# filtering benchmarks by topic, which we'll cover in a later section.
SUITE["utf8"] = BenchmarkGroup(["string", "unicode"])
SUITE["trig"] = BenchmarkGroup(["math", "triangles"])
SUITE["Substitution Models"] = BenchmarkGroup(["Trait", "Nucleic Acid"])

# Add some benchmarks to the "utf8" group
teststr = join(rand(MersenneTwister(1), 'a':'d', 10^4));
SUITE["utf8"]["replace"] = @benchmarkable replace($teststr, "a" => "b")
SUITE["utf8"]["join"] = @benchmarkable join($teststr, $teststr)

# Add some benchmarks to the "trig" group
for f in (sin, cos, tan)
    for x in (0.0, pi)
        SUITE["trig"][string(f), x] = @benchmarkable $(f)($x)
    end
end

#Substitution Models
m1 = EqualRatesSubstitutionModel(2, [1.0], ["low","high"]);
SUITE["Substitution Models"]["Trait"] = @benchmarkable P!(P(m1, 1.0), m1, 3.0)
m1 = JC69([0.25])
SUITE["Substitution Models"]["Nucleic Acid"] = @benchmarkable P!(P(m1, 1.0), m1, 3.0)

# #uploads to Github
# using GitHub, JSON, PkgBenchmark

# results = benchmarkpkg("PkgBenchmark");

# gist_json = JSON.parse(
#             """
#             {
#             "description": "A benchmark for PkgBenchmark",
#             "public": false,
#             "files": {
#                 "benchmark.md": {
#                 "content": "$(escape_string(sprint(export_markdown, results)))"
#                 }
#             }
#             }
#             """
#         )

# posted_gist = create_gist(params = gist_json);

# url = get(posted_gist.html_url)
# #URI(https://gist.github.com/317378b4fcf2fb4c5585b104c3b177a8) #TODO