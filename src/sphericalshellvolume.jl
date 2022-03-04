"""

```
sphericalshellvolume(i,step)
```

Computes the volume of the spherical shell defined within [(i-1)*step,i*step].

"""
function sphericalshellvolume(i, step)
    rmin = (i - 1) * step
    return (4 * pi / 3) * ((rmin + step)^3 - rmin^3)
end
