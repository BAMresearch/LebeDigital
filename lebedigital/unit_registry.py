import pint

ureg = pint.UnitRegistry(cache_folder=":auto:")   # initialize unit registry

# user defined dimensions
ureg.define('[moment] = [force] * [length]')
ureg.define('[stress] = [force] / [length]**2')
ureg.define('[line_load] = [force] / [length]'
