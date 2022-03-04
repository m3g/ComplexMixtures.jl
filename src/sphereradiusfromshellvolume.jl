"""

```
sphereradiusfromshellvolume(volume,step)
```

Computes the radius that corresponds to a spherical shell of
a given volume.

"""
function sphereradiusfromshellvolume(volume, step)

    fourthirdsofpi = 4 * pi / 3
    if 3 * step * volume - pi * step^4 <= 0.0
        return 0.0
    end
    rmin =
        (sqrt(3 * pi) * sqrt(3 * step * volume - pi * step^4) - 3 * pi * step^2) /
        (6 * pi * step)
    return (0.5 * (volume / fourthirdsofpi + 2 * rmin^3))^(1 / 3)

end
