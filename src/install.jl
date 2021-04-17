using Pkg

add(x) = Pkg.add(x)

function main()
    dependencies = [
        "Plots",
        "GaussianProcesses"
    ]
    for dep in dependencies
        add(dep)
    end
end

main()