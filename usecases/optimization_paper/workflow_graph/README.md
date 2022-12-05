## Design variable = [a,b .. ] with description of a,b
* Slag to concrete ratio : value between 0 and 1 describing the amount of slag compared to cement
* aggregate to water+cement ratio : value between 0 and 1 describing amount of aggregates
## Random Variable = [c,d .. ] with description
I guess everything based on the data from the parameter identification, mainly:
* paste compressive strength
* paste Youngs modulus
* hydration parameters

but from these follow
* concrete compressive strength
* concrete Youngs modulus
* concrete tensile strength
* some more homogenization output values
* fem outputs
* beam design output
## Objectives : and how it relates to the design variable and random variable (for example time(a,c))
* minimize GWP
less concrete should improve GWP, therefore
* more aggregates should reduce GWP
* more slag should reduce GWP
## Constrains: and how it relates to the design variable and random variable
* beam design: output needs to be positive (sufficient concrete strength for given steel and load)
  * more aggregates should reduce strength
  * more slag should reduce strength
* max temperature: max temperature should not exceed limit
  * more aggregates should reduce temperature output
  * more slag should reduce temperature output
* time of demolding: time when stress yield limit is not exceeded
  * difficult to judge, as complex relation
  * I guess more aggregates could increase the time as temperature is lowered and strength
  * I guess more slag could increase the time as temperature is lowered and strength reduced