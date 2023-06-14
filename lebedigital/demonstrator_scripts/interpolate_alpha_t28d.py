from lebedigital.unit_registry import ureg


@ureg.check("", "", "")
def interpolate_alpha_t28d(alpha_mix1, alpha_mix2, ratio):
    """
    This function linearly interpolates the alpha value at 28 days based on the slag ratio.
    In the future it might be worth it, to use a different interpolation function and include uncertainties.

    Parameters
    ----------
    ratio : float / pint unitless
        amount of slag compared to cement, value from 0 to 1
    alpha_mix1 : float / pint unitless
        degree of hydration after 28 days for mix 1 (OPC)
    alpha_mix2 : float / pint unitless
        degree of hydration after 28 days for mix 2 (100% slag)

    Returns
    -------
    alpha_t28d : float / pint unitless
        degree of hydration after 28 days for the mix with the given slag ratio
    """
    assert 0 <= ratio <= 1, "slag_ratio must be between 0 and 1"
    alpha_t28d = alpha_mix1 + (alpha_mix2 - alpha_mix1) * ratio

    return alpha_t28d
