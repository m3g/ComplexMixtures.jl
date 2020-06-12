
using DelimitedFiles
using Plots
nogtk()

old = readdlm("./gmd.dat",comments=true,comment_char='#')

plot(layout=(6,1))

sp=1
plot!(ylabel="MDDF",subplot=sp)
plot!(old[:,1],old[:,2],subplot=sp,label="old")
plot!(R.d,R.mddf,subplot=sp,label="new")

sp=2
plot!(ylabel="KB",subplot=sp)
plot!(old[:,1],old[:,3],subplot=sp,label="old")
plot!(R.d,R.kb,subplot=sp,label="new - rnd")
plot!(R.d,R.kb_shell,subplot=sp,label="new - shell")

sp=3
plot!(ylabel="Count",subplot=sp)
plot!(old[:,1],old[:,4],subplot=sp,label="old")
plot!(R.d,R.count,subplot=sp,label="new")
plot!(old[:,1],old[:,5],subplot=sp,label="old - rand")
plot!(R.d,R.count_random,subplot=sp,label="new -rand")

sp=4
plot!(ylabel="Shell vol", subplot=sp)
plot!(old[:,1],old[:,8],subplot=sp,label="old")
plot!(R.d,R.volume.shell,subplot=sp,label="new")

sp=5
plot!(ylabel="Sum", subplot=sp)
plot!(old[:,1],old[:,6],subplot=sp,label="old - md")
scatter!(R.d,R.sum_count,subplot=sp,label="new - md")

sp=6
plot!(old[:,1],old[:,7],subplot=sp,label="old - rand")
scatter!(R.d,R.sum_count_random,subplot=sp,label="new - rand")
scatter!(R.d,R.sum_shell,subplot=sp,label="new - shell")

plot!(size=(600,1600))
savefig("./plots.pdf")
