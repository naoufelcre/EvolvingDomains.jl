## skeleton of a typical problem we want to attack with the package ! 

### init

geometry = field -> a level set field ! 
state = field -> a scalar field representing state ! (technically could be any kind of field: vector valued, tensor valued)

-> DYNAMICS PROBLEM
function dynamics(config) # can take anything really the point is that it returns a value, it can depend on what you want it to depend on !
	return systemRate
end

### solving loop

for or while time condition unsatisfied

	currentRate = dynamics(config)
	currentVelocity = velFromRates(currentRate)

	geometry = eulerianSolverStep(currentVelocity) # can be done with trixi or package built-in
	state = LagrangianSolverStep(currentVelocity) # or not ? for instance a purely diffusive state do not need this ! 

	updateConfig()
	
end


## some specifics

### kinematics toolbox

### dynamics toolbox

### other

mainly multilevel, the idea is that you can use it to explore ! 

