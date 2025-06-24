using CSV, DataFrames, SummaryTables, Typst_jll

datadir = joinpath(@__DIR__, "../data/")
club_data_fn = joinpath(datadir, "Club_Data_Table.csv")
coords = joinpath(datadir, "coordinate_speeds_clubgolfswing.csv")

df = CSV.read(club_data_fn, DataFrame)
tbl = simple_table(df)
io = IOBuffer()

show(stdout, MIME"text/typst"(), tbl)


df = CSV.read(coords, DataFrame)