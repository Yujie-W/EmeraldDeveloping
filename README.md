# Emerald

Vegetation makes the land a green piece out of the blue marble (Earth). Therefore, we call the model Emerald, standing for: **E**arth **M**odeling with **E**nhanced **R**epresentations of **A**reo-Aqua and **L**and **D**ynamics.

## Difference with Emerald.jl

This repository is an under continuous development, and therefore may have unexpected results, particularly for the code not covered in the UI and CI tests. For a more stable release, we refer the users to [**Emerald.jl**](https://github.com/silicormosia/Emerald.jl.git), where we keep the stable versions used in formal publications. However, both models are named **Emerald** and would conflict if you include both in your project.

Therefore, if you need the most recent features, you are recommended to use EmeraldDeveloping (this repo). If you need a stable version for general production, [**Emerald.jl**](https://github.com/silicormosia/Emerald.jl.git) is a better option.

## Model components
- `Land` Land Surface Model based on CliMA Land (Remastered)
- `Ocean` Ocean Model (RT under development)

## Installation
`Emerald` is not a registered package, to use `Emerald`, one need to install it manually through Julia REPL:
```julia
using Pkg;
Pkg.add(PackageSpec(url="https://github.com/Yujie-W/EmeraldDeveloping.git"));
```

To use more testing features, you may install a specific branch, say `wyujie`
```julia
using Pkg;
Pkg.add(PackageSpec(url="https://github.com/Yujie-W/EmeraldDeveloping.git", rev="wyujie"));
```

## Tutorials
For the details about how to use Emerald, you may look into the `test/tutorial` folder, and subscribe [`Bilibili`](https://space.bilibili.com/22576885) and [`Bilibili Live`](http://live.bilibili.com/14243417).
