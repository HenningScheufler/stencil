## stencil improvements for of2106


* subsetStencil -> creates stencil on a subset of all cells e.g. interfaceregion
* stencilLooper -> 

```
    stencilLooper<scalar> looper = maskedStencil.stencilValues(cellNum);
    count = 0;
    forAll(looper,i)
    {
        for (const label& val:looper[i])
        {
            count += val;
        }
    }
```

* improved performance of cellToCellStencil by a factor of 6

This is relevant for  simulation with topological changes e.g. AMR

```
damBreakWithObstacle interIsoFoam

new:
    ExecutionTime = 164.29 s  ClockTime = 164 s

    End

    Finalising parallel run

old:
    ExecutionTime = 188.39 s  ClockTime = 189 s

    End

    Finalising parallel run

```