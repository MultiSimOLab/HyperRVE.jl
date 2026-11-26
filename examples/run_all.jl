# Runs every example in its own module; exits with an error if any fails.
failed = 0
for f in sort(filter(f -> endswith(f, ".jl") && basename(f) != "run_all.jl", readdir(@__DIR__; join=true)))
    print(rpad(basename(f), 36))
    t = @elapsed try
        m = Module()
        Core.eval(m, :(include(x) = Base.include($m, x)))
        redirect_stdout(devnull) do
            Core.eval(m, :(display(x) = nothing))
            Base.include(m, f)
        end
        print("OK   ")
    catch e
        global failed += 1
        print("FAIL ", sprint(showerror, e)[1:min(end, 200)])
    end
    println(round(t; digits=1), " s")
end
failed == 0 || error("$failed examples failed")
