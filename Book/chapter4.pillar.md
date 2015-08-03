

## Multi\-Pathogen/Multi\-Host Models



###1\.  SIR model with three species of hosts

In the standard models of epidemiology, the population is compartmentalized by only clinic status\.
As such, the population has only one degree of subdivision\.
In a context of multi\-host \(multi\-species\) model, the host population has two degrees of subdivision due to the attribute species of each individual\.



####1\.1\.  Equations




####1\.2\.  Configurations of Kendrick model

In epidemiology, it is important to distinguish between two basic assumptions in terms of the underlying structure of contacts within the population\.
Either the model is assumed to be mass action or pseudo mass action\.
The first kind reflects the situation where the number of contacts is independent of the population size\.
So that the force of infection \.
In some circumstances, the transmission rate  is rescaled by \.
The second one assumes that as the population size increases, so does the contact rate\.
As such the force of infection \.

At the moment, Kendrick model includes three parameters of configuration: sizeOfPopulation, rescale, mass\_action\.
By default:



```smalltalk
#sizeOfPopulation->#population
#rescale->true
#mass_action->true
```



In the context of the multi\-species model, it is important to config the size of population for each species\.
As such:



```smalltalk
model configurations: {#sizeOfPopulation->#(#species)}
```





####1\.3\.  Kendrick model

In this model, we define the parameter  for three scopes corresponding to each species\.
In order to represent the interaction between three species, we define a contact network\.
Due to this network, the force of infection will be modified as:

where  denotes the strength of connection between species  and \.



```smalltalk
| model graph |
	model := KEModel new.
	model
		population:
			(KEMetaPopulation new
				attributes:
					{(#status -> #(#S #I #R)).
					(#species -> #(#mosquito #reservoir1 #reservoir2))}).
	model
		buildFromAttributes: #(#status #species)
		compartments:
			{(#(#S #mosquito) -> 9800).
			(#(#I #mosquito) -> 200).
			(#(#R #mosquito) -> 0).
			(#(#S #reservoir1) -> 1000).
			(#(#I #reservoir1) -> 0).
			(#(#R #reservoir1) -> 0).
			(#(#S #reservoir2) -> 2000).
			(#(#I #reservoir2) -> 0).
			(#(#R #reservoir2) -> 0)}.
	model addParameter: #mu
		   inScopes: {
				#species->#mosquito.
				#species->#reservoir1.
				#species->#reservoir2}
		   values: #(12.17 0.05 0.05).
	model addParameter: #gamma value: 52.
	model addParameter: #beta value: 1.
	model addParameter: #N value: #sizeOfPopulation.
	model configurations: { #sizeOfPopulation->#(#species) }.

	graph := KEContactNetwork
			newOn: model population
			atAttribute: #species.
	graph edges: { #mosquito->#reservoir1. #mosquito->#reservoir2 };
			strengthOfAllConnections: 0.02.
	model
		addTransitionFrom: '{#status: #S}'
		to: '{#status: #I}'
		probability: [ :m | (m atParameter: #beta) * (m probabilityOfContact: '{#status: #I}') ].
	model
		addTransitionFrom: '{#status: #I}'
		to: '{#status: #R}'
		probability: [ :m | m atParameter: #gamma ].
	model
		addTransitionFrom: '{#status: #S}'
		to: #empty
		probability: [ :m | m atParameter: #mu ].
	model
		addTransitionFrom: '{#status: #I}'
		to: #empty
		probability: [ :m | m atParameter: #mu ].
	model
		addTransitionFrom: '{#status: #R}'
		to: #empty
		probability: [ :m | m atParameter: #mu ].
	model
		addTransitionFrom: #empty
		to: '{#status: #S}'
		probability: [ :m | m atParameter: #mu ].
```


