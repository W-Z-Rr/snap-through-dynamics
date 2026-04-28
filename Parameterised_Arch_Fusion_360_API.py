"""Fusion 360 script for generating an arched beam centreline."""

import math
import traceback

import adsk.core
import adsk.fusion


# Initialise the global Application and UserInterface objects.
app = adsk.core.Application.get()
ui = app.userInterface

def run(_context: str):
    """Called by Fusion 360 when the script is run."""

    try:
        plot_arch_shape_by_span()

    except:  # pylint: disable=bare-except
        app.log(f'Failed:\n{traceback.format_exc()}')

def plot_arch_shape_by_span():
    """Create a fitted spline for the arched beam profile on the XY plane."""

    ui_local = None

    try:
        app_local = adsk.core.Application.get()
        ui_local = app_local.userInterface

        design = adsk.fusion.Design.cast(app_local.activeProduct)
        root = design.rootComponent

        # Create a new sketch on the XY plane.
        sketch = root.sketches.add(root.xYConstructionPlane)

        # Arch inputs

        desired_span_B = 15          # Physical span between ends (cm)

        beta = 0.01                  # End angle parameter
        c_const = 2                  # Shape parameter

        num_points = 100             # Number of integration points

        # Non-dimensional arch shape

        nondim_points = []

        s_start = 0.0
        s_end = 1.0
        delta_s = (s_end - s_start) / num_points

        x_nondim = 0.0
        y_nondim = 0.0

        nondim_points.append((x_nondim,y_nondim))

        for i in range(num_points):
            s = s_start + i * delta_s

            term1 = (1 - 2*s) * beta
            term2 = c_const * s * (1 - s) * (1 - 2*s)
            theta = term1 + term2

            dx_nondim = math.cos(theta) * delta_s
            dy_nondim = math.sin(theta) * delta_s

            x_nondim += dx_nondim
            y_nondim += dy_nondim

            nondim_points.append((x_nondim,y_nondim))

        # Scale to physical span

        b_nondim = nondim_points[-1][0]

        scale_factor = 1.0

        if b_nondim != 0:
            scale_factor = desired_span_B / b_nondim
        else:
            ui_local.messageBox('Error: non-dimensional span is zero. Cannot scale.')
            return

        # Create sketch points

        points = adsk.core.ObjectCollection.create()

        for x_nd, y_nd in nondim_points:
            x_physical = x_nd * scale_factor
            y_physical = y_nd * scale_factor

            points.add(adsk.core.Point3D.create(x_physical,y_physical,0))

        # Create fitted spline

        sketch.sketchCurves.sketchFittedSplines.add(points)

    except:  # pylint: disable=bare-except
        if ui_local:
            ui_local.messageBox('Failed:\n{}'.format(traceback.format_exc()))

# Run the function directly when the script is executed.
plot_arch_shape_by_span()
