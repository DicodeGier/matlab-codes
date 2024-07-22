%{
Hoi Rens, ik hoop je met dit voorbeeld te laten zien hoe Gurobi werkt,
althans hetgene wat je vermoedelijk nodig gaat hebben. Als er iets niet
duidelijk is, app met gerust!

-Dico
%}

%{
Om het simpel te houden, laat ik het probleem klein. Als je de basis
begrijpt kun je het daarna zelf uitbreiden. Ik ga programmeren:

max x_1 + 2x_2
s.t.
    3x_1 + 5x_2 <= 20
    2x_1 <= 3x_2 - 5
    x_1 >= 0
    x_2 >= 0
    x_1, x_2 integer
%}

clear
clc

%De volgende regel kopieer ik er standaard in om gurobi te gebruiken (ik ben vergeten of dit nodig is), meer
%hierover in de meegestuurde manual
addpath("C:\gurobi1001\win64\matlab")

%{
Nu volgt het belangrijkste punt van Gurobi, je moet in je objective
function, constraints en nog een aantal dingen IEDERE BESLISVARIABELE
SPECIFICEREN. In ons probleem moet ik dus voor de objective function
en iedere constraint iets specificeren voor x_1 en x_2. Als je dus de
constraint 2x_1 <= 3x_2 - 5 wilt programmeren zul je hem dus moeten omschrijven naar 2x_1 - 3x_2 <=
-5 omdat je in de constraint matrix een waarde voor x_2 moet hebben, dit kun
je (voor zover ik weet) niet aan de rechterkant van het '<=' teken laten staan.
%}

%initialiseer het model
model = [];

%{
Aan dit model moeten de volgende dingen worden toegevoegd:
-obj: objective function
-A: constraint matrix
-rhs: right hand side van de constraints
-sense: de (on)gelijkheid per constraint
-modelsense: minimaliseren of maximaliseren
-vtype: type variable (integer, continuous, binary, ...)

je kunt ook nog andere dingen toevoegen zoals lower- en upperbounds per
variabele maar dat laat ik hier achterwege. Je kan dat uiteraard ook gewoon
als constraint toevoegen zoals ik ga doen
%}

%objective function
objective_function = [1 2]; %deze geef ik als rijvector mee
model.obj = objective_function; %voeg toe aan het model

%constraint matrix
constraint_matrix = [3 5;
                     2 -3; %constraint die we hebben omgeschreven
                     1 0; % x_1 + 0x_2 voor x_1 >= 0
                     0 1]; % 0x_1 + x_2 voor x_2 >= 0
%{
Gurobi accepteert bovenstaande matrix NIET, het moet namelijk een sparse
matrix zijn, een matrix waarin alleen de niet-0 elementen staan. Als je
matrix niet te groot is kun je dit simpelweg oplossen met
sparse(constraint_matrix), maar als je echt heel heel veel constraints of
beslisvariabelen hebt kan die memory issues krijgen, hier zou ik niet
vanuit gaan, ik heb heel ORM met sparse(...) kunnen doen
%}
sparse_constraint_matrix = sparse(constraint_matrix);
model.A = sparse_constraint_matrix;

%rhs
right_hand_side = [20;-5;0;0]; %deze geef ik als colomvector mee
model.rhs = right_hand_side;

%sense
sense = ['<';'<';'>';'>']; % de opties zijn <,> of = en ik geef hem als colomvector mee
%{
met maar 4 constraints kun je dit handmatig uittypen, In ons geval 
is bv. de command sense = [repmat('<',2,1) repmat('>',2,1]) handig waarbij
2 het aantal keer is dat de ongelijkheid achter elkaar door wordt gebruikt
%}
model.sense = sense;

%modelsense
%uiteraard kun je dingen ook meteen aan het model toevoegen zonder
%tussenvariabelen
model.modelsense = 'max';

%vtype
%degene die ik tot nu toe gebruikt heb zijn I: integer, C: continuous, B:
%binary
model.vtype = ['I','I']; %deze geef ik als rijvector mee
%uiteraard kan ook hier repmat handig zijn met bv 'I' ipv '<'

%{
Gurobi geeft (zelfs met een ; erachter) best veel output in de command
window, voor een keertje maakt dat niet uit maar als je je programma 200
keer gaat oproepen achter elkaar is dat toch niet zo heel prettig om te
zien en al helemaal omdat de uiteindelijke struct toch in je workspace
komt te staan, vandaar de volgende command (ik begrijp ook niet helemaal
wat het doet, credits naar Daan van Turnhout)
%}
params.outputflag = 0;
 

%{
Gurobi heeft nu alles wat die nodig heeft, dus je kan hem oplossen.
Normaal is alleen result = gurobi(model) genoeg, maar omdat we de output
willen onderdrukken moeten we een kleine toevoeging doen
%}

result = gurobi(model,params);

%{
De uitkomst staat nu in de struct genaamd result, met verschillende fields
waarbij vooral de objval en x voor ons van belang zijn, deze kun je ook heel
eenvoudig uit de struct halen
%}

optimal_value = result.objval; %want de fieldname in result is objval
optimal_choices = result.x; %want de fieldname is x

%hierdoor kan je het resultaat gelijk netjes printen

fprintf('the optimal value for x1 = %g, for x2 = %g and the optimal value is %g\n',optimal_choices(1),optimal_choices(2),optimal_value)
fprintf('ik weet ook niet waarom die nou -0 zegt ipv 0 maar maakt denk ik niet zo veel uit\n')
