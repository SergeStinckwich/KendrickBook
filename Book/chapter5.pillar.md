

## Spatial models



###1\.  SEIR model with spatial dynamics

We investigate the impact of spatial effects\.
Considering a spatial model organised by n patches arranged in a ring\.
The individuals can move between two neighbouring patches\.
In each patch, we have a sub\-population with four compartments S, E, I and R\.
To specify this model, we use the Migration Network built in Kendrick\. Due to this network, the model will have migration transitions from one compartment to other\.



####1\.1\.  Equations





####1\.2\.  Kendrick model



```smalltalk
| model graph |
model := KEModel new.
model population: KEMetaPopulation new.
model population attributes: {
	#patch->((1 to: 5)).
	#status->#(S E I R)}.
model
	buildFromAttributes: #(#status #patch)
	compartments: {
		  (#(#S 1) -> 900). (#(#E 1) -> 0). (#(#I 1) -> 100). (#(#R 1) -> 0).
        (#(#S 2) -> 1000). (#(#E 2) -> 0). (#(#I 2) -> 0). (#(#R 2) -> 0).
        (#(#S 3) -> 1000). (#(#E 3) -> 0). (#(#I 3) -> 0). (#(#R 3) -> 0).
        (#(#S 4) -> 1000). (#(#E 4) -> 0). (#(#I 4) -> 0). (#(#R 4) -> 0).
        (#(#S 5) -> 1000). (#(#E 5) -> 0). (#(#I 5) -> 0). (#(#R 5) -> 0).
	}.
model
	addParameter: #beta
	inScopes: {
		(#patch->1).
		(#patch->2).
		(#patch->3).
		(#patch->4).
		(#patch->5)
	}
	values: #(0.75 0.5 0.5 0.5 0.5).
model addParameter: '{
	#v: 0.00274,
	#d: 0.0000365,
	#epsilon: 0.5,
	#gamma: 0.25,
	#N: #sizeOfPopulation}'.
model configurations: {
		#sizeOfPopulation->#(#patch).
		#rescale->false }.
graph := KEMigrationNetwork
				newOn: model population
				atAttribute: #patch
				topology: (KERandomSmallWorldNetwork new K: 2; beta: 0).
graph strengthOfAllConnections: 0.03.
graph addMigrationTransitionsTo: model.

model addEquation: 'S:t=d*N-d*S-beta*S*I+v*R' parseAsAnEquation.
model addEquation: 'E:t=beta*S*I-d*E-epsilon*E' parseAsAnEquation.
model addEquation: 'I:t=epsilon*E-d*I-gamma*I' parseAsAnEquation.
model addEquation: 'R:t=gamma*I-d*R-v*R' parseAsAnEquation.
```


