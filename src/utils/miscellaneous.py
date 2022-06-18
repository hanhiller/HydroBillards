import typer
import numpy as np


def highlighted_str(value: str, prefix: str = "", suffix: str = ""):
    return " ".join([prefix, typer.style(value, fg=typer.colors.GREEN, bold=True), suffix])


def calc_emin(external_particle_density, simulating_particle_density, temp):
    # TODO: Understand these constants and put them in proper place
    return 1 - 78312.7 / ((np.sqrt(external_particle_density) * simulating_particle_density) / temp ** 1.5) ** 1.25
