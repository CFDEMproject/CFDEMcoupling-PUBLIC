particle_mesh nGridPoints 10

particle_data number_particles 1005 
#coupling json myCoupling
coupling liggghts

#Heat properties
model propertiesThermo heatThermalConductivity_solid 
model propertiesThermo heatCapacity_solid 
model propertiesThermo heatDensity_solid 


#Equations
#Option 1: cool with implicit coupling
#modelEqn 1DSpherical  heat     BC0 1  BC1 2	# 2 Convective: will cool based on coefficient and fluid temp

#Option 2: cool with explicit setting of flux
modelEqn 1DSpherical  heat     BC0 1  BC1 0	# 0 Neumann: will cool with a fixed rate provided by LIGGGHTS

control outputTimeStep 0.1
control timeStep 5e-6
#control run 0.000025 init yes
