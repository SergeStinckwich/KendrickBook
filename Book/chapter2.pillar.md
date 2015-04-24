

## Introduction to Simple Epidemic Model
<a name=" Introduction to Simple Epidemic Model"></a>


###1\.  Simple SIR \(without births and deaths\)
<a name=" Simple SIR (without births and deaths)"></a>
Program 2\.1 is a simple SIR model \(page 19 of the book\)\. These are the equations and the code of the model:



####1\.1\.  Equations
<a name=" Equations"></a>



####1\.2\.  Pharo code
<a name=" Pharo code"></a>

```smalltalk
|solver system dt beta gamma values stepper diag colors maxTime|
dt := 1.0.
beta := 1.4247.
gamma := 0.14286.
maxTime := 70.0.
system := ExplicitSystem block: [ :x :t| |c|
     c := Array new: 3.
     c at: 1 put: (beta negated) * (x at: 1) * (x at: 2).
     c at: 2 put: (beta * (x at: 1) * (x at: 2)) - (gamma * (x at: 2)).
     c at: 3 put: gamma * (x at: 2).
     c
     ].
stepper := RungeKuttaStepper onSystem: system.
solver := (ExplicitSolver new) stepper: stepper; system: system; dt: dt.
state := { 1-1e-6. 1e-6. 0}.
values := (0.0 to: maxTime by: dt) collect: [ :t| |state| state := stepper doStep: state
                                                          time: t stepSize: dt ].
diag := OrderedCollection new.
colors := Array with: Color blue with: Color red with: Color green.
1 to: 3 do: [ :i|
    diag add:
        ((GETLineDiagram new)
            models: (1 to: maxTime+1 by: 1);
            y: [ :x| (values at: x) at: i ];
            color: (colors at: i))
     ].
builder := (GETDiagramBuilder new).
builder compositeDiagram
    xAxisLabel: 'Time in days';
    yAxisLabel: 'Number of Individuals';
    regularAxis;
    diagrams: diag.
builder open.
```





###2\.  SIR model with births and deaths
<a name=" SIR model with births and deaths"></a>

####2\.1\.  Equations
<a name=" Equations"></a>



####2\.2\.  Pharo code
<a name=" Pharo code"></a>

```smalltalk
|solver system dt beta gamma values stepper diag mu colors maxTime|
dt := 1.0.
mu := 1/(70*365.0).
beta := 520/365.0.
gamma := 1/7.0.
maxTime := 60*365.
system := ExplicitSystem block: [ :x :t| |c|
     c := Array new: 3.
     c at: 1 put: mu - (beta  * (x at: 1) * (x at: 2)) - (mu * (x at:1)).
     c at: 2 put: (beta * (x at: 1) * (x at: 2)) - (gamma * (x at: 2)) - (mu * (x at:2)).
     c at: 3 put: (gamma * (x at: 2)) - (mu * (x at: 2)).
     c
     ].

stepper := RungeKuttaStepper onSystem: system.
solver := (ExplicitSolver new) stepper: stepper; system: system; dt: dt.
state := { 0.1. 1e-4. 1-0.1-1e-4}.
values := (0.0 to: maxTime by: dt) collect: [ :t| |state| state := stepper doStep: state
                                                          time: t stepSize: dt ].
diag := OrderedCollection new.
colors := Array with: Color blue with: Color red with: Color green.
1 to: 3 do: [ :i|
    diag add:
        ((GETLineDiagram new)
            models: (1 to: maxTime+1 by: 1);
            y: [ :x| (values at: x) at: i ];
            color: (colors at: i))
     ].
builder := (GETDiagramBuilder new).
builder compositeDiagram
    xAxisLabel: 'Time in days';
    yAxisLabel: 'Number of Individuals';
    regularAxis;
    diagrams: diag.
builder open.
```





###3\.  SIR model with disease induced mortality and density dependent transmission
<a name=" SIR model with disease induced mortality and density dependent transmission"></a>

####3\.1\.  Equations
<a name=" Equations"></a>



###4\.  SIR model, disease induced mortality and frequency dependent transmission
<a name=" SIR model, disease induced mortality and frequency dependent transmission"></a>

####4\.1\.  Equations
<a name=" Equations"></a>


####4\.2\.  Pharo code
<a name=" Pharo code"></a>


###5\.  SIS model without births or deaths
<a name=" SIS model without births or deaths"></a>

####5\.1\. Equations
<a name="Equations"></a>



###6\.  SEIR model with births and deaths
<a name=" SEIR model with births and deaths"></a>

###7\.  SIR with a carrier state
<a name=" SIR with a carrier state"></a>