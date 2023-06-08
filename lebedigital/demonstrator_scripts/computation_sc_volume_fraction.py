import pytest

from lebedigital.unit_registry import ureg


def computation_sc_volume_fraction(input_dic):
    # initialize output dictionary
    output = {}

    if input_dic["sc_mass_fraction"] == 0:
        return 0
    else:
        mass_sub = 1 * ureg("kg") / (1 + (1 - input_dic["sc_mass_fraction"]) / input_dic["sc_mass_fraction"])
        print(f"mass_sub = {mass_sub}")
        mass_cem = 1 * ureg("kg") - mass_sub
        print(f"mass_cem = {mass_cem}")
        vol_sub = mass_sub / input_dic["density_sub"]
        vol_cem = mass_cem / input_dic["density_cem"]
        sc_volume_fraction = vol_sub / (vol_sub + vol_cem)

        return sc_volume_fraction


# run when script is call
if __name__ == "__main__":
    # input dictionary
    input_dic = {
        "sc_mass_fraction": 0.7,
        "density_sub": 100 * ureg("kg/m^3"),
        "density_cem": 600 * ureg("kg/m^3"),
    }

    # compute the volume fraction of the substitute
    sc_volume_fraction = computation_sc_volume_fraction(input_dic)
    print(f"sc_volume_fraction = {sc_volume_fraction}")
