import spotpy
import numpy as np

# Dict([["α_cp", [0.4, 1.0]],
#                  ["μ1", [0.0, 0.4]],
#                  ["ω", [0.0, 1.0]],
#                  ["ϵ", [0.0, 0.1]],
#                  ["κ_upper", [0.0, 0.05]],
#                  ["ψ_E", [0.0, 0.1]],
#                  ["ψ_Q", [0.0, 0.5]],
#                  ["ψ_P", [0.0, 0.5]],
#                  ["p_f", [0.0, 1.0]]])


class spotpy_setup(object):
    def __init__(self):
        self.parameters = [
            spotpy.parameter.Uniform('μ1', 0.0, 0.4),
            spotpy.parameter.Uniform('ω', 0.0, 1.0),
            spotpy.parameter.Uniform('ϵ', 0.0, 0.1),
            spotpy.parameter.Uniform('κ_upper', [0.0, 0.05])
        ]

    def parameters(self):
        return spotpy.parameter.generate(self.parameters)

    # def 


if __name__ == "__main__":
    calibrate_parameters()