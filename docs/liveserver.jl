
import Pkg
Pkg.activate("LiveServer"; shared = true)
Pkg.add("LiveServer")

import LiveServer

Pkg.activate(@__DIR__)
LiveServer.servedocs()
