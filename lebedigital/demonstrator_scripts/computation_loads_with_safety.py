from lebedigital.unit_registry import ureg


@ureg.check("", "", "N/m", "N/m", "N", "N")
def computation_loads_with_safety(
    safety_factor_permanent,
    safety_factor_variable,
    distributed_load_permanent,
    distributed_load_variable,
    point_load_permanent,
    point_load_variable,
):
    """
    This function computes the loads with safety factors.

    Args:
        safety_factor_permanent:
        safety_factor_variable:
        distributed_load_permanent:
        distributed_load_variable:
        point_load_permanent:
        point_load_variable:

    Returns:
        distributed_load:
        point_load:
    """
    distributed_load = (
        distributed_load_permanent * safety_factor_permanent + distributed_load_variable * safety_factor_variable
    )
    point_load = point_load_permanent * safety_factor_permanent + point_load_variable * safety_factor_variable

    return distributed_load, point_load
