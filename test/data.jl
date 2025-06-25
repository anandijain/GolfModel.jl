using CSV, DataFrames, SummaryTables, Typst_jll, Typstry

datadir = joinpath(@__DIR__, "../data/")
club_data_fn = joinpath(datadir, "Club_Data_Table.csv")
coords = joinpath(datadir, "coordinate_speeds_clubgolfswing.csv")

df = CSV.read(club_data_fn, DataFrame)
tbl = simple_table(df)
io = IOBuffer()

show(stdout, MIME"text/typst"(), tbl)

df = CSV.read(coords, DataFrame)
df_orig = deepcopy(df)
sw_start = 3.667
sw_end = 3.967
sw_end - sw_start
df = filter(x -> x.time > 3.5, df)
df


getd(d, xs) = map(x -> d[x], xs)
sym_order = [t, th1, om1, th2, om2, l1, m1, l2, m2, g, tau_sh, tau_wr]
# making a table for a system 
function sys2df(sys)
    syms = model_syms(sys)
    nms = getdescription.(syms)
    f(x) = "`$x`"
    strsyms = f.(string.(syms))
    typst_syms = ["\$t\$", "\$theta_1\$", "\$omega_1\$",
        "\$theta_2\$", "\$omega_2\$",
        "\$l_1\$", "\$m_1\$", "\$l_2\$", "\$m_2\$",
        "\$g\$", "\$tau_(\"shoulder\")\$", "\$tau_(\"wrist\")\$"]
    units = ModelingToolkit.get_unit.(syms)
    defvals = getd(defs, syms)
    model_table = DataFrame(
        Symbol("Code Symbol") => strsyms,
        Symbol("Paper Symbol") => typst_syms,
        :Name => nms,
        :Unit => units,
        :Value => defvals,
    )

end
mdl_tbl = simple_table(model_table)
show(stdout, MIME"text/typst"(), mdl_tbl)