# Spline visualization
The splines are plotted using `gustaf` and `vedo` library.

The following will give a brief introduction into spline visualization and its options.

## Creating a basic NURBS
![NURBS](../source/_static/readme_nurbs.png)
For example here we will be creating nurbs with array-like inputs.
```python
import splinepy

# Initialize nurbs with any array-like input
nurbs = splinepy.NURBS(
    degrees=[2, 1],
    knot_vectors=[
        [0, 0, 0, 1, 1, 1],
        [0, 0, 1, 1],
    ],
    control_points=[
        [-1.0, 0.0],
        [-1.0, 1.0],
        [0.0, 1.0],
        [-2.0, 0.0],
        [-2.0, 2.0],
        [0.0, 2.0],
    ],
    weights=[
        [1.0],
        [2**-0.5],
        [1.0],
        [1.0],
        [2**-0.5],
        [1.0],
    ],
)

# vizusalize
nurbs.show()
```


## Setting show_options
![plotting_show_options](../source/_static/plotting_show_options.png)
You can also customize visualization by setting ***show_options***:
```python
nurbs.show_options["c"] = "pink"
nurbs.show_options["control_point_ids"] = False

# You can use `update()` to set multiple options at once
# As visualization runs through mesh based engines, you can set sampling
# resolutions per dimension.
# In case of int, it will apply that for all dimension
nurbs.show_options.update(control_mesh_c="green", resolutions=[201, 3])

nurbs.show()
```

## Plotting data on spline
![plotting_data](../source/_static/plotting_data.png)

This does not need to be only from 2D->3D but can also be from 1D->2D etc.

This is just the beginning of what you can do with the given plotting and creation capabilities in this library. Please look into the <a href='https://github.com/tataratat/splinepy/tree/main/examples'>examples folder</a> to see more capabilities.


## Notebook plotting

You can also plot your splines inside a notebook. For this you need to change the vedo backend to 'k3d'. To do this you need to add the following lines to the top of the notebook.

```
import vedo

vedo.settings.default_backend = "k3d"
```

After this most functionality of the normal plotting should be available to you. Please write an issue if you find something that you can plot in a script but not in a notebook. There is one known issue regarding number of Line objects that we can plot (limited to 200 from `vedo`). We are working on a solution, but meanwhile, plotting curves is limited to resolutions of 201; for 2D/3D splines you can turn off `knots` in `show_options`. An example is provided in the folder `examples/ipynb`.
