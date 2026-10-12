using Documenter
using Emerald


# define default docs pages
pages = Pair{Any,Any}[
    "Home" => "index.md",
];


# format the docs
mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(),
        )
    )
);

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1
);


# build the docs
makedocs(
    sitename = "EmeraldDeveloping",
    format = format,
    clean = false,
    modules = [Emerald],
    pages = pages,
    warnonly = [:cross_references, :missing_docs],
);


# deploy the docs to Github gh-pages
deploydocs(
    repo = "github.com/Yujie-W/EmeraldDeveloping.git",
    target = "build",
    devbranch = "main",
    push_preview = true
);
