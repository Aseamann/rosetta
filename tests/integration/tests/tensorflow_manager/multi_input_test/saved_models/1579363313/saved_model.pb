Õ±

õ%Į%
:
Add
x"T
y"T
z"T"
Ttype:
2	
P
Assert
	condition
	
data2T"
T
list(type)(0"
	summarizeint
E
AssignAddVariableOp
resource
value"dtype"
dtypetype
B
AssignVariableOp
resource
value"dtype"
dtypetype
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
~
BiasAddGrad
out_backprop"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
R
BroadcastGradientArgs
s0"T
s1"T
r0"T
r1"T"
Ttype0:
2	
N
Cast	
x"SrcT	
y"DstT"
SrcTtype"
DstTtype"
Truncatebool( 
I
ConcatOffset

concat_dim
shape*N
offset*N"
Nint(0
h
ConcatV2
values"T*N
axis"Tidx
output"T"
Nint(0"	
Ttype"
Tidxtype0:
2	
8
Const
output"dtype"
valuetensor"
dtypetype
¹
DenseToDenseSetOperation	
set1"T	
set2"T
result_indices	
result_values"T
result_shape	"
set_operationstring"
validate_indicesbool("
Ttype:
	2	
5
DivNoNan
x"T
y"T
z"T"
Ttype:
2
S
DynamicStitch
indices*N
data"T*N
merged"T"
Nint(0"	
Ttype
B
Equal
x"T
y"T
z
"
Ttype:
2	

W

ExpandDims

input"T
dim"Tdim
output"T"	
Ttype"
Tdimtype0:
2	
^
Fill
dims"
index_type

value"T
output"T"	
Ttype"

index_typetype0:
2	
?
FloorDiv
x"T
y"T
z"T"
Ttype:
2	
9
FloorMod
x"T
y"T
z"T"
Ttype:

2	
=
Greater
x"T
y"T
z
"
Ttype:
2	
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
;
Maximum
x"T
y"T
z"T"
Ttype:

2	

Mean

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
N
Merge
inputs"T*N
output"T
value_index"	
Ttype"
Nint(0
=
Mul
x"T
y"T
z"T"
Ttype:
2	
.
Neg
x"T
y"T"
Ttype:

2	

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
X
PlaceholderWithDefault
input"dtype
output"dtype"
dtypetype"
shapeshape
6
Pow
x"T
y"T
z"T"
Ttype:

2	

Prod

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
~
RandomUniform

shape"T
output"dtype"
seedint "
seed2int "
dtypetype:
2"
Ttype:
2	
a
Range
start"Tidx
limit"Tidx
delta"Tidx
output"Tidx"
Tidxtype0:	
2	
@
ReadVariableOp
resource
value"dtype"
dtypetype
>
RealDiv
x"T
y"T
z"T"
Ttype:
2	
E
Relu
features"T
activations"T"
Ttype:
2	
V
ReluGrad
	gradients"T
features"T
	backprops"T"
Ttype:
2	
[
Reshape
tensor"T
shape"Tshape
output"T"	
Ttype"
Tshapetype0:
2	
ą
ResourceApplyAdam
var
m
v
beta1_power"T
beta2_power"T
lr"T

beta1"T

beta2"T
epsilon"T	
grad"T" 
Ttype:
2	"
use_lockingbool( "
use_nesterovbool( 
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
P
Shape

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
e
ShapeN
input"T*N
output"out_type*N"
Nint(0"	
Ttype"
out_typetype0:
2	
0
Sigmoid
x"T
y"T"
Ttype:

2
=
SigmoidGrad
y"T
dy"T
z"T"
Ttype:

2
O
Size

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
a
Slice

input"T
begin"Index
size"Index
output"T"	
Ttype"
Indextype:
2	
G
SquaredDifference
x"T
y"T
z"T"
Ttype:

2	
:
Sub
x"T
y"T
z"T"
Ttype:
2	

Sum

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
M
Switch	
data"T
pred

output_false"T
output_true"T"	
Ttype
c
Tile

input"T
	multiples"
Tmultiples
output"T"	
Ttype"

Tmultiplestype0:
2	
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape
9
VarIsInitializedOp
resource
is_initialized
"train*2.0.0-alpha02v1.12.0-9492-g2c319fb4158ŲĆ
f
in1Placeholder*
dtype0*'
_output_shapes
:’’’’’’’’’*
shape:’’’’’’’’’
f
in2Placeholder*
dtype0*'
_output_shapes
:’’’’’’’’’*
shape:’’’’’’’’’
§
1in1_dense/kernel/Initializer/random_uniform/shapeConst*
valueB"      *#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes
:

/in1_dense/kernel/Initializer/random_uniform/minConst*
valueB
 *×³]æ*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes
: 

/in1_dense/kernel/Initializer/random_uniform/maxConst*
valueB
 *×³]?*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes
: 
Ų
9in1_dense/kernel/Initializer/random_uniform/RandomUniformRandomUniform1in1_dense/kernel/Initializer/random_uniform/shape*
T0*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes

:
Ž
/in1_dense/kernel/Initializer/random_uniform/subSub/in1_dense/kernel/Initializer/random_uniform/max/in1_dense/kernel/Initializer/random_uniform/min*
_output_shapes
: *
T0*#
_class
loc:@in1_dense/kernel
š
/in1_dense/kernel/Initializer/random_uniform/mulMul9in1_dense/kernel/Initializer/random_uniform/RandomUniform/in1_dense/kernel/Initializer/random_uniform/sub*
T0*#
_class
loc:@in1_dense/kernel*
_output_shapes

:
ā
+in1_dense/kernel/Initializer/random_uniformAdd/in1_dense/kernel/Initializer/random_uniform/mul/in1_dense/kernel/Initializer/random_uniform/min*
_output_shapes

:*
T0*#
_class
loc:@in1_dense/kernel
”
in1_dense/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape
:*!
shared_namein1_dense/kernel*#
_class
loc:@in1_dense/kernel
q
1in1_dense/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpin1_dense/kernel*
_output_shapes
: 

in1_dense/kernel/AssignAssignVariableOpin1_dense/kernel+in1_dense/kernel/Initializer/random_uniform*#
_class
loc:@in1_dense/kernel*
dtype0

$in1_dense/kernel/Read/ReadVariableOpReadVariableOpin1_dense/kernel*
dtype0*
_output_shapes

:*#
_class
loc:@in1_dense/kernel

 in1_dense/bias/Initializer/zerosConst*
valueB*    *!
_class
loc:@in1_dense/bias*
dtype0*
_output_shapes
:

in1_dense/biasVarHandleOp*
dtype0*
_output_shapes
: *
shape:*
shared_namein1_dense/bias*!
_class
loc:@in1_dense/bias
m
/in1_dense/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpin1_dense/bias*
_output_shapes
: 

in1_dense/bias/AssignAssignVariableOpin1_dense/bias in1_dense/bias/Initializer/zeros*!
_class
loc:@in1_dense/bias*
dtype0

"in1_dense/bias/Read/ReadVariableOpReadVariableOpin1_dense/bias*!
_class
loc:@in1_dense/bias*
dtype0*
_output_shapes
:
p
in1_dense/MatMul/ReadVariableOpReadVariableOpin1_dense/kernel*
dtype0*
_output_shapes

:
r
in1_dense/MatMulMatMulin1in1_dense/MatMul/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’
k
 in1_dense/BiasAdd/ReadVariableOpReadVariableOpin1_dense/bias*
dtype0*
_output_shapes
:

in1_dense/BiasAddBiasAddin1_dense/MatMul in1_dense/BiasAdd/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’
[
in1_dense/ReluReluin1_dense/BiasAdd*
T0*'
_output_shapes
:’’’’’’’’’
§
1in2_dense/kernel/Initializer/random_uniform/shapeConst*
valueB"   
   *#
_class
loc:@in2_dense/kernel*
dtype0*
_output_shapes
:

/in2_dense/kernel/Initializer/random_uniform/minConst*
dtype0*
_output_shapes
: *
valueB
 *ó5æ*#
_class
loc:@in2_dense/kernel

/in2_dense/kernel/Initializer/random_uniform/maxConst*
dtype0*
_output_shapes
: *
valueB
 *ó5?*#
_class
loc:@in2_dense/kernel
Ų
9in2_dense/kernel/Initializer/random_uniform/RandomUniformRandomUniform1in2_dense/kernel/Initializer/random_uniform/shape*
T0*#
_class
loc:@in2_dense/kernel*
dtype0*
_output_shapes

:

Ž
/in2_dense/kernel/Initializer/random_uniform/subSub/in2_dense/kernel/Initializer/random_uniform/max/in2_dense/kernel/Initializer/random_uniform/min*
T0*#
_class
loc:@in2_dense/kernel*
_output_shapes
: 
š
/in2_dense/kernel/Initializer/random_uniform/mulMul9in2_dense/kernel/Initializer/random_uniform/RandomUniform/in2_dense/kernel/Initializer/random_uniform/sub*
T0*#
_class
loc:@in2_dense/kernel*
_output_shapes

:

ā
+in2_dense/kernel/Initializer/random_uniformAdd/in2_dense/kernel/Initializer/random_uniform/mul/in2_dense/kernel/Initializer/random_uniform/min*
T0*#
_class
loc:@in2_dense/kernel*
_output_shapes

:

”
in2_dense/kernelVarHandleOp*!
shared_namein2_dense/kernel*#
_class
loc:@in2_dense/kernel*
dtype0*
_output_shapes
: *
shape
:

q
1in2_dense/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpin2_dense/kernel*
_output_shapes
: 

in2_dense/kernel/AssignAssignVariableOpin2_dense/kernel+in2_dense/kernel/Initializer/random_uniform*#
_class
loc:@in2_dense/kernel*
dtype0

$in2_dense/kernel/Read/ReadVariableOpReadVariableOpin2_dense/kernel*
dtype0*
_output_shapes

:
*#
_class
loc:@in2_dense/kernel

 in2_dense/bias/Initializer/zerosConst*
valueB
*    *!
_class
loc:@in2_dense/bias*
dtype0*
_output_shapes
:


in2_dense/biasVarHandleOp*!
_class
loc:@in2_dense/bias*
dtype0*
_output_shapes
: *
shape:
*
shared_namein2_dense/bias
m
/in2_dense/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpin2_dense/bias*
_output_shapes
: 

in2_dense/bias/AssignAssignVariableOpin2_dense/bias in2_dense/bias/Initializer/zeros*
dtype0*!
_class
loc:@in2_dense/bias

"in2_dense/bias/Read/ReadVariableOpReadVariableOpin2_dense/bias*!
_class
loc:@in2_dense/bias*
dtype0*
_output_shapes
:

p
in2_dense/MatMul/ReadVariableOpReadVariableOpin2_dense/kernel*
dtype0*
_output_shapes

:

r
in2_dense/MatMulMatMulin2in2_dense/MatMul/ReadVariableOp*'
_output_shapes
:’’’’’’’’’
*
T0
k
 in2_dense/BiasAdd/ReadVariableOpReadVariableOpin2_dense/bias*
dtype0*
_output_shapes
:


in2_dense/BiasAddBiasAddin2_dense/MatMul in2_dense/BiasAdd/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’

[
in2_dense/ReluReluin2_dense/BiasAdd*
T0*'
_output_shapes
:’’’’’’’’’

S
merge/concat/axisConst*
value	B :*
dtype0*
_output_shapes
: 

merge/concatConcatV2in1_dense/Reluin2_dense/Relumerge/concat/axis*
T0*
N*'
_output_shapes
:’’’’’’’’’
”
.dense3/kernel/Initializer/random_uniform/shapeConst*
dtype0*
_output_shapes
:*
valueB"      * 
_class
loc:@dense3/kernel

,dense3/kernel/Initializer/random_uniform/minConst*
dtype0*
_output_shapes
: *
valueB
 *qÄæ* 
_class
loc:@dense3/kernel

,dense3/kernel/Initializer/random_uniform/maxConst*
valueB
 *qÄ?* 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes
: 
Ļ
6dense3/kernel/Initializer/random_uniform/RandomUniformRandomUniform.dense3/kernel/Initializer/random_uniform/shape*
T0* 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes

:
Ņ
,dense3/kernel/Initializer/random_uniform/subSub,dense3/kernel/Initializer/random_uniform/max,dense3/kernel/Initializer/random_uniform/min*
T0* 
_class
loc:@dense3/kernel*
_output_shapes
: 
ä
,dense3/kernel/Initializer/random_uniform/mulMul6dense3/kernel/Initializer/random_uniform/RandomUniform,dense3/kernel/Initializer/random_uniform/sub*
T0* 
_class
loc:@dense3/kernel*
_output_shapes

:
Ö
(dense3/kernel/Initializer/random_uniformAdd,dense3/kernel/Initializer/random_uniform/mul,dense3/kernel/Initializer/random_uniform/min*
T0* 
_class
loc:@dense3/kernel*
_output_shapes

:

dense3/kernelVarHandleOp* 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes
: *
shape
:*
shared_namedense3/kernel
k
.dense3/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense3/kernel*
_output_shapes
: 

dense3/kernel/AssignAssignVariableOpdense3/kernel(dense3/kernel/Initializer/random_uniform* 
_class
loc:@dense3/kernel*
dtype0

!dense3/kernel/Read/ReadVariableOpReadVariableOpdense3/kernel* 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes

:

dense3/bias/Initializer/zerosConst*
valueB*    *
_class
loc:@dense3/bias*
dtype0*
_output_shapes
:

dense3/biasVarHandleOp*
_class
loc:@dense3/bias*
dtype0*
_output_shapes
: *
shape:*
shared_namedense3/bias
g
,dense3/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense3/bias*
_output_shapes
: 

dense3/bias/AssignAssignVariableOpdense3/biasdense3/bias/Initializer/zeros*
_class
loc:@dense3/bias*
dtype0

dense3/bias/Read/ReadVariableOpReadVariableOpdense3/bias*
_class
loc:@dense3/bias*
dtype0*
_output_shapes
:
j
dense3/MatMul/ReadVariableOpReadVariableOpdense3/kernel*
dtype0*
_output_shapes

:
u
dense3/MatMulMatMulmerge/concatdense3/MatMul/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’
e
dense3/BiasAdd/ReadVariableOpReadVariableOpdense3/bias*
dtype0*
_output_shapes
:
y
dense3/BiasAddBiasAdddense3/MatMuldense3/BiasAdd/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’
[
dense3/SigmoidSigmoiddense3/BiasAdd*'
_output_shapes
:’’’’’’’’’*
T0

dense3_targetPlaceholder*%
shape:’’’’’’’’’’’’’’’’’’*
dtype0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’
R
ConstConst*
valueB*  ?*
dtype0*
_output_shapes
:

dense3_sample_weightsPlaceholderWithDefaultConst*
dtype0*#
_output_shapes
:’’’’’’’’’*
shape:’’’’’’’’’
v
total/Initializer/zerosConst*
valueB
 *    *
_class

loc:@total*
dtype0*
_output_shapes
: 
x
totalVarHandleOp*
dtype0*
_output_shapes
: *
shape: *
shared_nametotal*
_class

loc:@total
[
&total/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal*
_output_shapes
: 
g
total/AssignAssignVariableOptotaltotal/Initializer/zeros*
_class

loc:@total*
dtype0
q
total/Read/ReadVariableOpReadVariableOptotal*
_class

loc:@total*
dtype0*
_output_shapes
: 
v
count/Initializer/zerosConst*
dtype0*
_output_shapes
: *
valueB
 *    *
_class

loc:@count
x
countVarHandleOp*
shape: *
shared_namecount*
_class

loc:@count*
dtype0*
_output_shapes
: 
[
&count/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount*
_output_shapes
: 
g
count/AssignAssignVariableOpcountcount/Initializer/zeros*
dtype0*
_class

loc:@count
q
count/Read/ReadVariableOpReadVariableOpcount*
_class

loc:@count*
dtype0*
_output_shapes
: 
\
metrics/accuracy/Cast/xConst*
valueB
 *   ?*
dtype0*
_output_shapes
: 
~
metrics/accuracy/GreaterGreaterdense3/Sigmoidmetrics/accuracy/Cast/x*
T0*'
_output_shapes
:’’’’’’’’’
z
metrics/accuracy/Cast_1Castmetrics/accuracy/Greater*

SrcT0
*'
_output_shapes
:’’’’’’’’’*

DstT0

metrics/accuracy/EqualEqualdense3_targetmetrics/accuracy/Cast_1*
T0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’

metrics/accuracy/Cast_2Castmetrics/accuracy/Equal*

SrcT0
*0
_output_shapes
:’’’’’’’’’’’’’’’’’’*

DstT0
r
'metrics/accuracy/Mean/reduction_indicesConst*
dtype0*
_output_shapes
: *
valueB :
’’’’’’’’’

metrics/accuracy/MeanMeanmetrics/accuracy/Cast_2'metrics/accuracy/Mean/reduction_indices*
T0*#
_output_shapes
:’’’’’’’’’
`
metrics/accuracy/ConstConst*
valueB: *
dtype0*
_output_shapes
:
k
metrics/accuracy/SumSummetrics/accuracy/Meanmetrics/accuracy/Const*
T0*
_output_shapes
: 
e
$metrics/accuracy/AssignAddVariableOpAssignAddVariableOptotalmetrics/accuracy/Sum*
dtype0

metrics/accuracy/ReadVariableOpReadVariableOptotal%^metrics/accuracy/AssignAddVariableOp^metrics/accuracy/Sum*
dtype0*
_output_shapes
: 
U
metrics/accuracy/SizeSizemetrics/accuracy/Mean*
T0*
_output_shapes
: 
f
metrics/accuracy/Cast_3Castmetrics/accuracy/Size*

SrcT0*
_output_shapes
: *

DstT0

&metrics/accuracy/AssignAddVariableOp_1AssignAddVariableOpcountmetrics/accuracy/Cast_3%^metrics/accuracy/AssignAddVariableOp*
dtype0
Æ
!metrics/accuracy/ReadVariableOp_1ReadVariableOpcount%^metrics/accuracy/AssignAddVariableOp'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

*metrics/accuracy/div_no_nan/ReadVariableOpReadVariableOptotal'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

,metrics/accuracy/div_no_nan/ReadVariableOp_1ReadVariableOpcount'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
¢
metrics/accuracy/div_no_nanDivNoNan*metrics/accuracy/div_no_nan/ReadVariableOp,metrics/accuracy/div_no_nan/ReadVariableOp_1*
T0*
_output_shapes
: 
c
metrics/accuracy/IdentityIdentitymetrics/accuracy/div_no_nan*
T0*
_output_shapes
: 
^
metrics/accuracy/Cast_4/xConst*
valueB
 *   ?*
dtype0*
_output_shapes
: 

metrics/accuracy/Greater_1Greaterdense3/Sigmoidmetrics/accuracy/Cast_4/x*
T0*'
_output_shapes
:’’’’’’’’’
|
metrics/accuracy/Cast_5Castmetrics/accuracy/Greater_1*

SrcT0
*'
_output_shapes
:’’’’’’’’’*

DstT0

metrics/accuracy/Equal_1Equaldense3_targetmetrics/accuracy/Cast_5*
T0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’

metrics/accuracy/Cast_6Castmetrics/accuracy/Equal_1*0
_output_shapes
:’’’’’’’’’’’’’’’’’’*

DstT0*

SrcT0

t
)metrics/accuracy/Mean_1/reduction_indicesConst*
valueB :
’’’’’’’’’*
dtype0*
_output_shapes
: 

metrics/accuracy/Mean_1Meanmetrics/accuracy/Cast_6)metrics/accuracy/Mean_1/reduction_indices*
T0*#
_output_shapes
:’’’’’’’’’
b
metrics/accuracy/Const_1Const*
valueB: *
dtype0*
_output_shapes
:
s
metrics/accuracy/Mean_2Meanmetrics/accuracy/Mean_1metrics/accuracy/Const_1*
T0*
_output_shapes
: 
¤
5loss/dense3_loss/mean_squared_error/SquaredDifferenceSquaredDifferencedense3/Sigmoiddense3_target*
T0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’

:loss/dense3_loss/mean_squared_error/Mean/reduction_indicesConst*
valueB :
’’’’’’’’’*
dtype0*
_output_shapes
: 
Ń
(loss/dense3_loss/mean_squared_error/MeanMean5loss/dense3_loss/mean_squared_error/SquaredDifference:loss/dense3_loss/mean_squared_error/Mean/reduction_indices*#
_output_shapes
:’’’’’’’’’*
T0
«
floss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapeShapedense3_sample_weights*
T0*
_output_shapes
:
§
eloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rankConst*
value	B :*
dtype0*
_output_shapes
: 
½
eloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapeShape(loss/dense3_loss/mean_squared_error/Mean*
T0*
_output_shapes
:
¦
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rankConst*
value	B :*
dtype0*
_output_shapes
: 
¦
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar/xConst*
value	B : *
dtype0*
_output_shapes
: 
Ł
bloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarEqualdloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar/xeloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank*
T0*
_output_shapes
: 
ć
nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/SwitchSwitchbloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarbloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar*
T0
*
_output_shapes
: : 

ploss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_tIdentityploss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch:1*
T0
*
_output_shapes
: 

ploss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_fIdentitynloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch*
T0
*
_output_shapes
: 

oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_idIdentitybloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar*
T0
*
_output_shapes
: 
é
ploss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1Switchbloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalaroloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
_output_shapes
: : *
T0
*u
_classk
igloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar
ė
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankEqualloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switchloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1*
_output_shapes
: *
T0

loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/SwitchSwitchdloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rankoloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
_output_shapes
: : *
T0*w
_classm
kiloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank

loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1Switcheloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rankoloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
T0*x
_classn
ljloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank*
_output_shapes
: : 
Ų
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/SwitchSwitchloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
T0
*
_output_shapes
: : 
Å
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_tIdentityloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:1*
T0
*
_output_shapes
: 
Ć
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_fIdentityloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch*
T0
*
_output_shapes
: 
Č
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_idIdentityloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
_output_shapes
: *
T0

ū
”loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dimConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
dtype0*
_output_shapes
: *
valueB :
’’’’’’’’’
¤
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims
ExpandDimsØloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1”loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim*
_output_shapes

:*
T0
¬
¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/SwitchSwitcheloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapeoloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id* 
_output_shapes
::*
T0*x
_classn
ljloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape

¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1Switch¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switchloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0*x
_classn
ljloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::

¢loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/ShapeConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
valueB"      *
dtype0*
_output_shapes
:
ó
¢loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/ConstConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
value	B :*
dtype0*
_output_shapes
: 

loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_likeFill¢loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Shape¢loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const*
T0*
_output_shapes

:
ļ
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axisConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
value	B :*
dtype0*
_output_shapes
: 
ø
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concatConcatV2loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDimsloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_likeloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis*
T0*
N*
_output_shapes

:
ż
£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dimConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
valueB :
’’’’’’’’’*
dtype0*
_output_shapes
: 
Ŗ
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1
ExpandDimsŖloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim*
T0*
_output_shapes

:
°
¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/SwitchSwitchfloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapeoloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
T0*y
_classo
mkloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::

Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1Switch¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switchloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0*y
_classo
mkloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::
å
«loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperationDenseToDenseSetOperationloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat*<
_output_shapes*
(:’’’’’’’’’:’’’’’’’’’:*
set_operationa-b*
T0
ż
£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dimsSize­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1*
_output_shapes
: *
T0
å
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/xConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
value	B : *
dtype0*
_output_shapes
: 
ś
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dimsEqualloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims*
T0*
_output_shapes
: 
ü
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1Switchloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0
*¤
_class
loc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
_output_shapes
: : 
ß
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/MergeMergeloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims*
T0
*
N*
_output_shapes
: : 
 
mloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/MergeMergeloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Mergerloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1*
T0
*
N*
_output_shapes
: : 
Ę
^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/ConstConst*8
value/B- B'weights can not be broadcast to values.*
dtype0*
_output_shapes
: 
Æ
`loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_1Const*
valueB Bweights.shape=*
dtype0*
_output_shapes
: 
ø
`loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_2Const*
dtype0*
_output_shapes
: *(
valueB Bdense3_sample_weights:0
®
`loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_3Const*
valueB Bvalues.shape=*
dtype0*
_output_shapes
: 
Ė
`loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_4Const*
dtype0*
_output_shapes
: *;
value2B0 B*loss/dense3_loss/mean_squared_error/Mean:0
«
`loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_5Const*
valueB B
is_scalar=*
dtype0*
_output_shapes
: 
ö
kloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/SwitchSwitchmloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Mergemloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge*
T0
*
_output_shapes
: : 

mloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_tIdentitymloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Switch:1*
T0
*
_output_shapes
: 

mloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_fIdentitykloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Switch*
T0
*
_output_shapes
: 

lloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_idIdentitymloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge*
T0
*
_output_shapes
: 
į
iloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/NoOpNoOpn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t

wloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependencyIdentitymloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_tj^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/NoOp*
_output_shapes
: *
T0
*
_classv
trloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t
Ź
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
dtype0*
_output_shapes
: *8
value/B- B'weights can not be broadcast to values.
±
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
valueB Bweights.shape=*
dtype0*
_output_shapes
: 
ŗ
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*(
valueB Bdense3_sample_weights:0*
dtype0*
_output_shapes
: 
°
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
valueB Bvalues.shape=*
dtype0*
_output_shapes
: 
Ķ
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*;
value2B0 B*loss/dense3_loss/mean_squared_error/Mean:0*
dtype0*
_output_shapes
: 
­
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
dtype0*
_output_shapes
: *
valueB B
is_scalar=


kloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/AssertAssertrloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switchrloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3*
T
2	

’
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/SwitchSwitchmloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Mergelloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
_output_shapes
: : *
T0
*
_classv
trloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge
ś
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1Switchfloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapelloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0*y
_classo
mkloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::
ų
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2Switcheloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapelloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0*x
_classn
ljloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::
ź
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3Switchbloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarlloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0
*u
_classk
igloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar*
_output_shapes
: : 

yloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1Identitymloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_fl^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert*
T0
*
_classv
trloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
_output_shapes
: 

jloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/MergeMergeyloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1wloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency*
T0
*
N*
_output_shapes
: : 

Sloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ShapeShape(loss/dense3_loss/mean_squared_error/Meank^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Merge*
T0*
_output_shapes
:

Sloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ConstConstk^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Merge*
valueB
 *  ?*
dtype0*
_output_shapes
: 
­
Mloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_likeFillSloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ShapeSloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/Const*
T0*#
_output_shapes
:’’’’’’’’’
Ž
Closs/dense3_loss/mean_squared_error/weighted_loss/broadcast_weightsMuldense3_sample_weightsMloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like*
T0*#
_output_shapes
:’’’’’’’’’
Ł
5loss/dense3_loss/mean_squared_error/weighted_loss/MulMul(loss/dense3_loss/mean_squared_error/MeanCloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights*#
_output_shapes
:’’’’’’’’’*
T0
`
loss/dense3_loss/ConstConst*
valueB: *
dtype0*
_output_shapes
:

loss/dense3_loss/SumSum5loss/dense3_loss/mean_squared_error/weighted_loss/Mulloss/dense3_loss/Const*
_output_shapes
: *
T0
}
loss/dense3_loss/num_elementsSize5loss/dense3_loss/mean_squared_error/weighted_loss/Mul*
T0*
_output_shapes
: 
y
"loss/dense3_loss/num_elements/CastCastloss/dense3_loss/num_elements*

SrcT0*
_output_shapes
: *

DstT0
[
loss/dense3_loss/mul/xConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 
x
loss/dense3_loss/mulMulloss/dense3_loss/mul/x"loss/dense3_loss/num_elements/Cast*
T0*
_output_shapes
: 
[
loss/dense3_loss/Const_1Const*
valueB *
dtype0*
_output_shapes
: 
n
loss/dense3_loss/Sum_1Sumloss/dense3_loss/Sumloss/dense3_loss/Const_1*
_output_shapes
: *
T0
q
loss/dense3_loss/valueDivNoNanloss/dense3_loss/Sum_1loss/dense3_loss/mul*
_output_shapes
: *
T0
O

loss/mul/xConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 
T
loss/mulMul
loss/mul/xloss/dense3_loss/value*
T0*
_output_shapes
: 
q
iter/Initializer/zerosConst*
dtype0	*
_output_shapes
: *
value	B	 R *
_class
	loc:@iter

iterVarHandleOp"/device:CPU:0*
dtype0	*
_output_shapes
: *
shape: *
shared_nameiter*
_class
	loc:@iter
h
%iter/IsInitialized/VarIsInitializedOpVarIsInitializedOpiter"/device:CPU:0*
_output_shapes
: 
r
iter/AssignAssignVariableOpiteriter/Initializer/zeros"/device:CPU:0*
_class
	loc:@iter*
dtype0	
}
iter/Read/ReadVariableOpReadVariableOpiter"/device:CPU:0*
_class
	loc:@iter*
dtype0	*
_output_shapes
: 

'learning_rate/Initializer/initial_valueConst*
valueB
 *o:* 
_class
loc:@learning_rate*
dtype0*
_output_shapes
: 

learning_rateVarHandleOp*
dtype0*
_output_shapes
: *
shape: *
shared_namelearning_rate* 
_class
loc:@learning_rate
k
.learning_rate/IsInitialized/VarIsInitializedOpVarIsInitializedOplearning_rate*
_output_shapes
: 

learning_rate/AssignAssignVariableOplearning_rate'learning_rate/Initializer/initial_value* 
_class
loc:@learning_rate*
dtype0

!learning_rate/Read/ReadVariableOpReadVariableOplearning_rate* 
_class
loc:@learning_rate*
dtype0*
_output_shapes
: 
~
decay/Initializer/initial_valueConst*
valueB
 *    *
_class

loc:@decay*
dtype0*
_output_shapes
: 
x
decayVarHandleOp*
shape: *
shared_namedecay*
_class

loc:@decay*
dtype0*
_output_shapes
: 
[
&decay/IsInitialized/VarIsInitializedOpVarIsInitializedOpdecay*
_output_shapes
: 
o
decay/AssignAssignVariableOpdecaydecay/Initializer/initial_value*
_class

loc:@decay*
dtype0
q
decay/Read/ReadVariableOpReadVariableOpdecay*
_class

loc:@decay*
dtype0*
_output_shapes
: 

 beta_1/Initializer/initial_valueConst*
valueB
 *fff?*
_class
loc:@beta_1*
dtype0*
_output_shapes
: 
{
beta_1VarHandleOp*
shared_namebeta_1*
_class
loc:@beta_1*
dtype0*
_output_shapes
: *
shape: 
]
'beta_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpbeta_1*
_output_shapes
: 
s
beta_1/AssignAssignVariableOpbeta_1 beta_1/Initializer/initial_value*
dtype0*
_class
loc:@beta_1
t
beta_1/Read/ReadVariableOpReadVariableOpbeta_1*
_class
loc:@beta_1*
dtype0*
_output_shapes
: 

 beta_2/Initializer/initial_valueConst*
dtype0*
_output_shapes
: *
valueB
 *w¾?*
_class
loc:@beta_2
{
beta_2VarHandleOp*
_class
loc:@beta_2*
dtype0*
_output_shapes
: *
shape: *
shared_namebeta_2
]
'beta_2/IsInitialized/VarIsInitializedOpVarIsInitializedOpbeta_2*
_output_shapes
: 
s
beta_2/AssignAssignVariableOpbeta_2 beta_2/Initializer/initial_value*
_class
loc:@beta_2*
dtype0
t
beta_2/Read/ReadVariableOpReadVariableOpbeta_2*
_class
loc:@beta_2*
dtype0*
_output_shapes
: 

!epsilon/Initializer/initial_valueConst*
valueB
 *æÖ3*
_class
loc:@epsilon*
dtype0*
_output_shapes
: 
~
epsilonVarHandleOp*
shape: *
shared_name	epsilon*
_class
loc:@epsilon*
dtype0*
_output_shapes
: 
_
(epsilon/IsInitialized/VarIsInitializedOpVarIsInitializedOpepsilon*
_output_shapes
: 
w
epsilon/AssignAssignVariableOpepsilon!epsilon/Initializer/initial_value*
dtype0*
_class
loc:@epsilon
w
epsilon/Read/ReadVariableOpReadVariableOpepsilon*
_class
loc:@epsilon*
dtype0*
_output_shapes
: 
`
training/Adam/gradients/ShapeConst*
valueB *
dtype0*
_output_shapes
: 
f
!training/Adam/gradients/grad_ys_0Const*
dtype0*
_output_shapes
: *
valueB
 *  ?

training/Adam/gradients/FillFilltraining/Adam/gradients/Shape!training/Adam/gradients/grad_ys_0*
T0*
_output_shapes
: 

)training/Adam/gradients/loss/mul_grad/MulMultraining/Adam/gradients/Fillloss/dense3_loss/value*
T0*
_output_shapes
: 
}
+training/Adam/gradients/loss/mul_grad/Mul_1Multraining/Adam/gradients/Fill
loss/mul/x*
_output_shapes
: *
T0
|
9training/Adam/gradients/loss/dense3_loss/value_grad/ShapeConst*
valueB *
dtype0*
_output_shapes
: 
~
;training/Adam/gradients/loss/dense3_loss/value_grad/Shape_1Const*
dtype0*
_output_shapes
: *
valueB 

Itraining/Adam/gradients/loss/dense3_loss/value_grad/BroadcastGradientArgsBroadcastGradientArgs9training/Adam/gradients/loss/dense3_loss/value_grad/Shape;training/Adam/gradients/loss/dense3_loss/value_grad/Shape_1*2
_output_shapes 
:’’’’’’’’’:’’’’’’’’’
®
>training/Adam/gradients/loss/dense3_loss/value_grad/div_no_nanDivNoNan+training/Adam/gradients/loss/mul_grad/Mul_1loss/dense3_loss/mul*
T0*
_output_shapes
: 
ź
7training/Adam/gradients/loss/dense3_loss/value_grad/SumSum>training/Adam/gradients/loss/dense3_loss/value_grad/div_no_nanItraining/Adam/gradients/loss/dense3_loss/value_grad/BroadcastGradientArgs*
T0*
_output_shapes
: 
Ū
;training/Adam/gradients/loss/dense3_loss/value_grad/ReshapeReshape7training/Adam/gradients/loss/dense3_loss/value_grad/Sum9training/Adam/gradients/loss/dense3_loss/value_grad/Shape*
_output_shapes
: *
T0
w
7training/Adam/gradients/loss/dense3_loss/value_grad/NegNegloss/dense3_loss/Sum_1*
_output_shapes
: *
T0
¼
@training/Adam/gradients/loss/dense3_loss/value_grad/div_no_nan_1DivNoNan7training/Adam/gradients/loss/dense3_loss/value_grad/Negloss/dense3_loss/mul*
T0*
_output_shapes
: 
Å
@training/Adam/gradients/loss/dense3_loss/value_grad/div_no_nan_2DivNoNan@training/Adam/gradients/loss/dense3_loss/value_grad/div_no_nan_1loss/dense3_loss/mul*
T0*
_output_shapes
: 
Ī
7training/Adam/gradients/loss/dense3_loss/value_grad/mulMul+training/Adam/gradients/loss/mul_grad/Mul_1@training/Adam/gradients/loss/dense3_loss/value_grad/div_no_nan_2*
T0*
_output_shapes
: 
ē
9training/Adam/gradients/loss/dense3_loss/value_grad/Sum_1Sum7training/Adam/gradients/loss/dense3_loss/value_grad/mulKtraining/Adam/gradients/loss/dense3_loss/value_grad/BroadcastGradientArgs:1*
_output_shapes
: *
T0
į
=training/Adam/gradients/loss/dense3_loss/value_grad/Reshape_1Reshape9training/Adam/gradients/loss/dense3_loss/value_grad/Sum_1;training/Adam/gradients/loss/dense3_loss/value_grad/Shape_1*
T0*
_output_shapes
: 

Atraining/Adam/gradients/loss/dense3_loss/Sum_1_grad/Reshape/shapeConst*
valueB *
dtype0*
_output_shapes
: 
ē
;training/Adam/gradients/loss/dense3_loss/Sum_1_grad/ReshapeReshape;training/Adam/gradients/loss/dense3_loss/value_grad/ReshapeAtraining/Adam/gradients/loss/dense3_loss/Sum_1_grad/Reshape/shape*
T0*
_output_shapes
: 
|
9training/Adam/gradients/loss/dense3_loss/Sum_1_grad/ConstConst*
valueB *
dtype0*
_output_shapes
: 
Ł
8training/Adam/gradients/loss/dense3_loss/Sum_1_grad/TileTile;training/Adam/gradients/loss/dense3_loss/Sum_1_grad/Reshape9training/Adam/gradients/loss/dense3_loss/Sum_1_grad/Const*
T0*
_output_shapes
: 

?training/Adam/gradients/loss/dense3_loss/Sum_grad/Reshape/shapeConst*
valueB:*
dtype0*
_output_shapes
:
ä
9training/Adam/gradients/loss/dense3_loss/Sum_grad/ReshapeReshape8training/Adam/gradients/loss/dense3_loss/Sum_1_grad/Tile?training/Adam/gradients/loss/dense3_loss/Sum_grad/Reshape/shape*
_output_shapes
:*
T0

7training/Adam/gradients/loss/dense3_loss/Sum_grad/ShapeShape5loss/dense3_loss/mean_squared_error/weighted_loss/Mul*
T0*
_output_shapes
:
ą
6training/Adam/gradients/loss/dense3_loss/Sum_grad/TileTile9training/Adam/gradients/loss/dense3_loss/Sum_grad/Reshape7training/Adam/gradients/loss/dense3_loss/Sum_grad/Shape*
T0*#
_output_shapes
:’’’’’’’’’
°
Xtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/ShapeShape(loss/dense3_loss/mean_squared_error/Mean*
T0*
_output_shapes
:
Ķ
Ztraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/Shape_1ShapeCloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights*
T0*
_output_shapes
:
ė
htraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/BroadcastGradientArgsBroadcastGradientArgsXtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/ShapeZtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/Shape_1*2
_output_shapes 
:’’’’’’’’’:’’’’’’’’’

Vtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/MulMul6training/Adam/gradients/loss/dense3_loss/Sum_grad/TileCloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights*
T0*#
_output_shapes
:’’’’’’’’’
Ā
Vtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/SumSumVtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/Mulhtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/BroadcastGradientArgs*
T0*
_output_shapes
:
Å
Ztraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/ReshapeReshapeVtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/SumXtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/Shape*
T0*#
_output_shapes
:’’’’’’’’’
ļ
Xtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/Mul_1Mul(loss/dense3_loss/mean_squared_error/Mean6training/Adam/gradients/loss/dense3_loss/Sum_grad/Tile*
T0*#
_output_shapes
:’’’’’’’’’
Č
Xtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/Sum_1SumXtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/Mul_1jtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/BroadcastGradientArgs:1*
T0*
_output_shapes
:
Ė
\training/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/Reshape_1ReshapeXtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/Sum_1Ztraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/Shape_1*
T0*#
_output_shapes
:’’’’’’’’’
°
Ktraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/ShapeShape5loss/dense3_loss/mean_squared_error/SquaredDifference*
T0*
_output_shapes
:
ģ
Jtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/SizeConst*
dtype0*
_output_shapes
: *
value	B :*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape
Ł
Itraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/addAdd:loss/dense3_loss/mean_squared_error/Mean/reduction_indicesJtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Size*
_output_shapes
: *
T0*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape
ķ
Itraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/modFloorModItraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/addJtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Size*
T0*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape*
_output_shapes
: 
š
Mtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape_1Const*
valueB *^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape*
dtype0*
_output_shapes
: 
ó
Qtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/range/startConst*
dtype0*
_output_shapes
: *
value	B : *^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape
ó
Qtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/range/deltaConst*
dtype0*
_output_shapes
: *
value	B :*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape
Ā
Ktraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/rangeRangeQtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/range/startJtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/SizeQtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/range/delta*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape*
_output_shapes
:
ņ
Ptraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Fill/valueConst*
value	B :*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape*
dtype0*
_output_shapes
: 
ō
Jtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/FillFillMtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape_1Ptraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Fill/value*
_output_shapes
: *
T0*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape
£
Straining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/DynamicStitchDynamicStitchKtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/rangeItraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/modKtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/ShapeJtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Fill*
N*
_output_shapes
:*
T0*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape
ń
Otraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Maximum/yConst*
dtype0*
_output_shapes
: *
value	B :*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape

Mtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/MaximumMaximumStraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/DynamicStitchOtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Maximum/y*
T0*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape*
_output_shapes
:
ū
Ntraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/floordivFloorDivKtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/ShapeMtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Maximum*
_output_shapes
:*
T0*^
_classT
RPloc:@training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape
Ä
Mtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/ReshapeReshapeZtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/weighted_loss/Mul_grad/ReshapeStraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/DynamicStitch*0
_output_shapes
:’’’’’’’’’’’’’’’’’’*
T0
¬
Jtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/TileTileMtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/ReshapeNtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/floordiv*
T0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’
²
Mtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape_2Shape5loss/dense3_loss/mean_squared_error/SquaredDifference*
T0*
_output_shapes
:
„
Mtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape_3Shape(loss/dense3_loss/mean_squared_error/Mean*
T0*
_output_shapes
:

Ktraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/ConstConst*
valueB: *
dtype0*
_output_shapes
:

Jtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/ProdProdMtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape_2Ktraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Const*
T0*
_output_shapes
: 

Mtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Const_1Const*
valueB: *
dtype0*
_output_shapes
:

Ltraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Prod_1ProdMtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Shape_3Mtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Const_1*
T0*
_output_shapes
: 

Qtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Maximum_1/yConst*
value	B :*
dtype0*
_output_shapes
: 

Otraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Maximum_1MaximumLtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Prod_1Qtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Maximum_1/y*
T0*
_output_shapes
: 

Ptraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/floordiv_1FloorDivJtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/ProdOtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Maximum_1*
T0*
_output_shapes
: 
Ō
Jtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/CastCastPtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/floordiv_1*

SrcT0*
_output_shapes
: *

DstT0
«
Mtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/truedivRealDivJtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/TileJtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/Cast*0
_output_shapes
:’’’’’’’’’’’’’’’’’’*
T0

Xtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/ShapeShapedense3/Sigmoid*
_output_shapes
:*
T0

Ztraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/Shape_1Shapedense3_target*
T0*
_output_shapes
:
ė
htraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/BroadcastGradientArgsBroadcastGradientArgsXtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/ShapeZtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/Shape_1*2
_output_shapes 
:’’’’’’’’’:’’’’’’’’’
ī
Ytraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/scalarConstN^training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/truediv*
valueB
 *   @*
dtype0*
_output_shapes
: 
Ā
Vtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/MulMulYtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/scalarMtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/truediv*
T0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’

Vtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/subSubdense3/Sigmoiddense3_targetN^training/Adam/gradients/loss/dense3_loss/mean_squared_error/Mean_grad/truediv*
T0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’
Ź
Xtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/mul_1MulVtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/MulVtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/sub*0
_output_shapes
:’’’’’’’’’’’’’’’’’’*
T0
Ä
Vtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/SumSumXtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/mul_1htraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/BroadcastGradientArgs*
T0*
_output_shapes
:
É
Ztraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/ReshapeReshapeVtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/SumXtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/Shape*'
_output_shapes
:’’’’’’’’’*
T0
Č
Xtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/Sum_1SumXtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/mul_1jtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/BroadcastGradientArgs:1*
T0*
_output_shapes
:
Ų
\training/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/Reshape_1ReshapeXtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/Sum_1Ztraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/Shape_1*
T0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’
ö
Vtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/NegNeg\training/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/Reshape_1*
T0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’
ä
7training/Adam/gradients/dense3/Sigmoid_grad/SigmoidGradSigmoidGraddense3/SigmoidZtraining/Adam/gradients/loss/dense3_loss/mean_squared_error/SquaredDifference_grad/Reshape*
T0*'
_output_shapes
:’’’’’’’’’
¤
7training/Adam/gradients/dense3/BiasAdd_grad/BiasAddGradBiasAddGrad7training/Adam/gradients/dense3/Sigmoid_grad/SigmoidGrad*
T0*
_output_shapes
:
×
1training/Adam/gradients/dense3/MatMul_grad/MatMulMatMul7training/Adam/gradients/dense3/Sigmoid_grad/SigmoidGraddense3/MatMul/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’*
transpose_b(
Ą
3training/Adam/gradients/dense3/MatMul_grad/MatMul_1MatMulmerge/concat7training/Adam/gradients/dense3/Sigmoid_grad/SigmoidGrad*
T0*
_output_shapes

:*
transpose_a(
p
.training/Adam/gradients/merge/concat_grad/RankConst*
value	B :*
dtype0*
_output_shapes
: 

-training/Adam/gradients/merge/concat_grad/modFloorModmerge/concat/axis.training/Adam/gradients/merge/concat_grad/Rank*
T0*
_output_shapes
: 
m
/training/Adam/gradients/merge/concat_grad/ShapeShapein1_dense/Relu*
T0*
_output_shapes
:

0training/Adam/gradients/merge/concat_grad/ShapeNShapeNin1_dense/Reluin2_dense/Relu*
T0*
N* 
_output_shapes
::

6training/Adam/gradients/merge/concat_grad/ConcatOffsetConcatOffset-training/Adam/gradients/merge/concat_grad/mod0training/Adam/gradients/merge/concat_grad/ShapeN2training/Adam/gradients/merge/concat_grad/ShapeN:1*
N* 
_output_shapes
::

/training/Adam/gradients/merge/concat_grad/SliceSlice1training/Adam/gradients/dense3/MatMul_grad/MatMul6training/Adam/gradients/merge/concat_grad/ConcatOffset0training/Adam/gradients/merge/concat_grad/ShapeN*
Index0*
T0*'
_output_shapes
:’’’’’’’’’

1training/Adam/gradients/merge/concat_grad/Slice_1Slice1training/Adam/gradients/dense3/MatMul_grad/MatMul8training/Adam/gradients/merge/concat_grad/ConcatOffset:12training/Adam/gradients/merge/concat_grad/ShapeN:1*
Index0*
T0*'
_output_shapes
:’’’’’’’’’

³
4training/Adam/gradients/in1_dense/Relu_grad/ReluGradReluGrad/training/Adam/gradients/merge/concat_grad/Slicein1_dense/Relu*'
_output_shapes
:’’’’’’’’’*
T0
µ
4training/Adam/gradients/in2_dense/Relu_grad/ReluGradReluGrad1training/Adam/gradients/merge/concat_grad/Slice_1in2_dense/Relu*
T0*'
_output_shapes
:’’’’’’’’’

¤
:training/Adam/gradients/in1_dense/BiasAdd_grad/BiasAddGradBiasAddGrad4training/Adam/gradients/in1_dense/Relu_grad/ReluGrad*
T0*
_output_shapes
:
¤
:training/Adam/gradients/in2_dense/BiasAdd_grad/BiasAddGradBiasAddGrad4training/Adam/gradients/in2_dense/Relu_grad/ReluGrad*
T0*
_output_shapes
:

Ś
4training/Adam/gradients/in1_dense/MatMul_grad/MatMulMatMul4training/Adam/gradients/in1_dense/Relu_grad/ReluGradin1_dense/MatMul/ReadVariableOp*'
_output_shapes
:’’’’’’’’’*
transpose_b(*
T0
·
6training/Adam/gradients/in1_dense/MatMul_grad/MatMul_1MatMulin14training/Adam/gradients/in1_dense/Relu_grad/ReluGrad*
T0*
_output_shapes

:*
transpose_a(
Ś
4training/Adam/gradients/in2_dense/MatMul_grad/MatMulMatMul4training/Adam/gradients/in2_dense/Relu_grad/ReluGradin2_dense/MatMul/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’*
transpose_b(
·
6training/Adam/gradients/in2_dense/MatMul_grad/MatMul_1MatMulin24training/Adam/gradients/in2_dense/Relu_grad/ReluGrad*
T0*
_output_shapes

:
*
transpose_a(
¼
2training/Adam/in1_dense/kernel/m/Initializer/zerosConst*
valueB*    *3
_class)
'%loc:@training/Adam/in1_dense/kernel/m*
dtype0*
_output_shapes

:
Ń
 training/Adam/in1_dense/kernel/mVarHandleOp*
dtype0*
_output_shapes
: *
shape
:*1
shared_name" training/Adam/in1_dense/kernel/m*3
_class)
'%loc:@training/Adam/in1_dense/kernel/m

Atraining/Adam/in1_dense/kernel/m/IsInitialized/VarIsInitializedOpVarIsInitializedOp training/Adam/in1_dense/kernel/m*
_output_shapes
: 
Ó
'training/Adam/in1_dense/kernel/m/AssignAssignVariableOp training/Adam/in1_dense/kernel/m2training/Adam/in1_dense/kernel/m/Initializer/zeros*3
_class)
'%loc:@training/Adam/in1_dense/kernel/m*
dtype0
Ź
4training/Adam/in1_dense/kernel/m/Read/ReadVariableOpReadVariableOp training/Adam/in1_dense/kernel/m*
dtype0*
_output_shapes

:*3
_class)
'%loc:@training/Adam/in1_dense/kernel/m
°
0training/Adam/in1_dense/bias/m/Initializer/zerosConst*
valueB*    *1
_class'
%#loc:@training/Adam/in1_dense/bias/m*
dtype0*
_output_shapes
:
Ē
training/Adam/in1_dense/bias/mVarHandleOp*/
shared_name training/Adam/in1_dense/bias/m*1
_class'
%#loc:@training/Adam/in1_dense/bias/m*
dtype0*
_output_shapes
: *
shape:

?training/Adam/in1_dense/bias/m/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/in1_dense/bias/m*
_output_shapes
: 
Ė
%training/Adam/in1_dense/bias/m/AssignAssignVariableOptraining/Adam/in1_dense/bias/m0training/Adam/in1_dense/bias/m/Initializer/zeros*1
_class'
%#loc:@training/Adam/in1_dense/bias/m*
dtype0
Ą
2training/Adam/in1_dense/bias/m/Read/ReadVariableOpReadVariableOptraining/Adam/in1_dense/bias/m*1
_class'
%#loc:@training/Adam/in1_dense/bias/m*
dtype0*
_output_shapes
:
¼
2training/Adam/in2_dense/kernel/m/Initializer/zerosConst*
valueB
*    *3
_class)
'%loc:@training/Adam/in2_dense/kernel/m*
dtype0*
_output_shapes

:

Ń
 training/Adam/in2_dense/kernel/mVarHandleOp*
shape
:
*1
shared_name" training/Adam/in2_dense/kernel/m*3
_class)
'%loc:@training/Adam/in2_dense/kernel/m*
dtype0*
_output_shapes
: 

Atraining/Adam/in2_dense/kernel/m/IsInitialized/VarIsInitializedOpVarIsInitializedOp training/Adam/in2_dense/kernel/m*
_output_shapes
: 
Ó
'training/Adam/in2_dense/kernel/m/AssignAssignVariableOp training/Adam/in2_dense/kernel/m2training/Adam/in2_dense/kernel/m/Initializer/zeros*3
_class)
'%loc:@training/Adam/in2_dense/kernel/m*
dtype0
Ź
4training/Adam/in2_dense/kernel/m/Read/ReadVariableOpReadVariableOp training/Adam/in2_dense/kernel/m*
dtype0*
_output_shapes

:
*3
_class)
'%loc:@training/Adam/in2_dense/kernel/m
°
0training/Adam/in2_dense/bias/m/Initializer/zerosConst*
dtype0*
_output_shapes
:
*
valueB
*    *1
_class'
%#loc:@training/Adam/in2_dense/bias/m
Ē
training/Adam/in2_dense/bias/mVarHandleOp*/
shared_name training/Adam/in2_dense/bias/m*1
_class'
%#loc:@training/Adam/in2_dense/bias/m*
dtype0*
_output_shapes
: *
shape:


?training/Adam/in2_dense/bias/m/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/in2_dense/bias/m*
_output_shapes
: 
Ė
%training/Adam/in2_dense/bias/m/AssignAssignVariableOptraining/Adam/in2_dense/bias/m0training/Adam/in2_dense/bias/m/Initializer/zeros*
dtype0*1
_class'
%#loc:@training/Adam/in2_dense/bias/m
Ą
2training/Adam/in2_dense/bias/m/Read/ReadVariableOpReadVariableOptraining/Adam/in2_dense/bias/m*1
_class'
%#loc:@training/Adam/in2_dense/bias/m*
dtype0*
_output_shapes
:

¶
/training/Adam/dense3/kernel/m/Initializer/zerosConst*
valueB*    *0
_class&
$"loc:@training/Adam/dense3/kernel/m*
dtype0*
_output_shapes

:
Č
training/Adam/dense3/kernel/mVarHandleOp*
shape
:*.
shared_nametraining/Adam/dense3/kernel/m*0
_class&
$"loc:@training/Adam/dense3/kernel/m*
dtype0*
_output_shapes
: 

>training/Adam/dense3/kernel/m/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense3/kernel/m*
_output_shapes
: 
Ē
$training/Adam/dense3/kernel/m/AssignAssignVariableOptraining/Adam/dense3/kernel/m/training/Adam/dense3/kernel/m/Initializer/zeros*
dtype0*0
_class&
$"loc:@training/Adam/dense3/kernel/m
Į
1training/Adam/dense3/kernel/m/Read/ReadVariableOpReadVariableOptraining/Adam/dense3/kernel/m*0
_class&
$"loc:@training/Adam/dense3/kernel/m*
dtype0*
_output_shapes

:
Ŗ
-training/Adam/dense3/bias/m/Initializer/zerosConst*
dtype0*
_output_shapes
:*
valueB*    *.
_class$
" loc:@training/Adam/dense3/bias/m
¾
training/Adam/dense3/bias/mVarHandleOp*.
_class$
" loc:@training/Adam/dense3/bias/m*
dtype0*
_output_shapes
: *
shape:*,
shared_nametraining/Adam/dense3/bias/m

<training/Adam/dense3/bias/m/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense3/bias/m*
_output_shapes
: 
æ
"training/Adam/dense3/bias/m/AssignAssignVariableOptraining/Adam/dense3/bias/m-training/Adam/dense3/bias/m/Initializer/zeros*.
_class$
" loc:@training/Adam/dense3/bias/m*
dtype0
·
/training/Adam/dense3/bias/m/Read/ReadVariableOpReadVariableOptraining/Adam/dense3/bias/m*.
_class$
" loc:@training/Adam/dense3/bias/m*
dtype0*
_output_shapes
:
¼
2training/Adam/in1_dense/kernel/v/Initializer/zerosConst*
dtype0*
_output_shapes

:*
valueB*    *3
_class)
'%loc:@training/Adam/in1_dense/kernel/v
Ń
 training/Adam/in1_dense/kernel/vVarHandleOp*
dtype0*
_output_shapes
: *
shape
:*1
shared_name" training/Adam/in1_dense/kernel/v*3
_class)
'%loc:@training/Adam/in1_dense/kernel/v

Atraining/Adam/in1_dense/kernel/v/IsInitialized/VarIsInitializedOpVarIsInitializedOp training/Adam/in1_dense/kernel/v*
_output_shapes
: 
Ó
'training/Adam/in1_dense/kernel/v/AssignAssignVariableOp training/Adam/in1_dense/kernel/v2training/Adam/in1_dense/kernel/v/Initializer/zeros*
dtype0*3
_class)
'%loc:@training/Adam/in1_dense/kernel/v
Ź
4training/Adam/in1_dense/kernel/v/Read/ReadVariableOpReadVariableOp training/Adam/in1_dense/kernel/v*3
_class)
'%loc:@training/Adam/in1_dense/kernel/v*
dtype0*
_output_shapes

:
°
0training/Adam/in1_dense/bias/v/Initializer/zerosConst*
valueB*    *1
_class'
%#loc:@training/Adam/in1_dense/bias/v*
dtype0*
_output_shapes
:
Ē
training/Adam/in1_dense/bias/vVarHandleOp*/
shared_name training/Adam/in1_dense/bias/v*1
_class'
%#loc:@training/Adam/in1_dense/bias/v*
dtype0*
_output_shapes
: *
shape:

?training/Adam/in1_dense/bias/v/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/in1_dense/bias/v*
_output_shapes
: 
Ė
%training/Adam/in1_dense/bias/v/AssignAssignVariableOptraining/Adam/in1_dense/bias/v0training/Adam/in1_dense/bias/v/Initializer/zeros*
dtype0*1
_class'
%#loc:@training/Adam/in1_dense/bias/v
Ą
2training/Adam/in1_dense/bias/v/Read/ReadVariableOpReadVariableOptraining/Adam/in1_dense/bias/v*1
_class'
%#loc:@training/Adam/in1_dense/bias/v*
dtype0*
_output_shapes
:
¼
2training/Adam/in2_dense/kernel/v/Initializer/zerosConst*
dtype0*
_output_shapes

:
*
valueB
*    *3
_class)
'%loc:@training/Adam/in2_dense/kernel/v
Ń
 training/Adam/in2_dense/kernel/vVarHandleOp*3
_class)
'%loc:@training/Adam/in2_dense/kernel/v*
dtype0*
_output_shapes
: *
shape
:
*1
shared_name" training/Adam/in2_dense/kernel/v

Atraining/Adam/in2_dense/kernel/v/IsInitialized/VarIsInitializedOpVarIsInitializedOp training/Adam/in2_dense/kernel/v*
_output_shapes
: 
Ó
'training/Adam/in2_dense/kernel/v/AssignAssignVariableOp training/Adam/in2_dense/kernel/v2training/Adam/in2_dense/kernel/v/Initializer/zeros*3
_class)
'%loc:@training/Adam/in2_dense/kernel/v*
dtype0
Ź
4training/Adam/in2_dense/kernel/v/Read/ReadVariableOpReadVariableOp training/Adam/in2_dense/kernel/v*3
_class)
'%loc:@training/Adam/in2_dense/kernel/v*
dtype0*
_output_shapes

:

°
0training/Adam/in2_dense/bias/v/Initializer/zerosConst*
dtype0*
_output_shapes
:
*
valueB
*    *1
_class'
%#loc:@training/Adam/in2_dense/bias/v
Ē
training/Adam/in2_dense/bias/vVarHandleOp*/
shared_name training/Adam/in2_dense/bias/v*1
_class'
%#loc:@training/Adam/in2_dense/bias/v*
dtype0*
_output_shapes
: *
shape:


?training/Adam/in2_dense/bias/v/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/in2_dense/bias/v*
_output_shapes
: 
Ė
%training/Adam/in2_dense/bias/v/AssignAssignVariableOptraining/Adam/in2_dense/bias/v0training/Adam/in2_dense/bias/v/Initializer/zeros*1
_class'
%#loc:@training/Adam/in2_dense/bias/v*
dtype0
Ą
2training/Adam/in2_dense/bias/v/Read/ReadVariableOpReadVariableOptraining/Adam/in2_dense/bias/v*
dtype0*
_output_shapes
:
*1
_class'
%#loc:@training/Adam/in2_dense/bias/v
¶
/training/Adam/dense3/kernel/v/Initializer/zerosConst*
dtype0*
_output_shapes

:*
valueB*    *0
_class&
$"loc:@training/Adam/dense3/kernel/v
Č
training/Adam/dense3/kernel/vVarHandleOp*
shape
:*.
shared_nametraining/Adam/dense3/kernel/v*0
_class&
$"loc:@training/Adam/dense3/kernel/v*
dtype0*
_output_shapes
: 

>training/Adam/dense3/kernel/v/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense3/kernel/v*
_output_shapes
: 
Ē
$training/Adam/dense3/kernel/v/AssignAssignVariableOptraining/Adam/dense3/kernel/v/training/Adam/dense3/kernel/v/Initializer/zeros*
dtype0*0
_class&
$"loc:@training/Adam/dense3/kernel/v
Į
1training/Adam/dense3/kernel/v/Read/ReadVariableOpReadVariableOptraining/Adam/dense3/kernel/v*0
_class&
$"loc:@training/Adam/dense3/kernel/v*
dtype0*
_output_shapes

:
Ŗ
-training/Adam/dense3/bias/v/Initializer/zerosConst*
valueB*    *.
_class$
" loc:@training/Adam/dense3/bias/v*
dtype0*
_output_shapes
:
¾
training/Adam/dense3/bias/vVarHandleOp*
shape:*,
shared_nametraining/Adam/dense3/bias/v*.
_class$
" loc:@training/Adam/dense3/bias/v*
dtype0*
_output_shapes
: 

<training/Adam/dense3/bias/v/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense3/bias/v*
_output_shapes
: 
æ
"training/Adam/dense3/bias/v/AssignAssignVariableOptraining/Adam/dense3/bias/v-training/Adam/dense3/bias/v/Initializer/zeros*.
_class$
" loc:@training/Adam/dense3/bias/v*
dtype0
·
/training/Adam/dense3/bias/v/Read/ReadVariableOpReadVariableOptraining/Adam/dense3/bias/v*
dtype0*
_output_shapes
:*.
_class$
" loc:@training/Adam/dense3/bias/v

9training/Adam/Adam/update_in1_dense/kernel/ReadVariableOpReadVariableOpiter"/device:CPU:0*
dtype0	*
_output_shapes
: 

0training/Adam/Adam/update_in1_dense/kernel/add/yConst*
dtype0	*
_output_shapes
: *
value	B	 R*#
_class
loc:@in1_dense/kernel
č
.training/Adam/Adam/update_in1_dense/kernel/addAdd9training/Adam/Adam/update_in1_dense/kernel/ReadVariableOp0training/Adam/Adam/update_in1_dense/kernel/add/y*
T0	*#
_class
loc:@in1_dense/kernel*
_output_shapes
: 
¼
/training/Adam/Adam/update_in1_dense/kernel/CastCast.training/Adam/Adam/update_in1_dense/kernel/add*

SrcT0	*#
_class
loc:@in1_dense/kernel*
_output_shapes
: *

DstT0
|
=training/Adam/Adam/update_in1_dense/kernel/Pow/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: 
ė
.training/Adam/Adam/update_in1_dense/kernel/PowPow=training/Adam/Adam/update_in1_dense/kernel/Pow/ReadVariableOp/training/Adam/Adam/update_in1_dense/kernel/Cast*
T0*#
_class
loc:@in1_dense/kernel*
_output_shapes
: 
~
?training/Adam/Adam/update_in1_dense/kernel/Pow_1/ReadVariableOpReadVariableOpbeta_2*
dtype0*
_output_shapes
: 
ļ
0training/Adam/Adam/update_in1_dense/kernel/Pow_1Pow?training/Adam/Adam/update_in1_dense/kernel/Pow_1/ReadVariableOp/training/Adam/Adam/update_in1_dense/kernel/Cast*
T0*#
_class
loc:@in1_dense/kernel*
_output_shapes
: 

Ktraining/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdam/ReadVariableOpReadVariableOplearning_rate*
dtype0*
_output_shapes
: 

Mtraining/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdam/ReadVariableOp_1ReadVariableOpbeta_1*
dtype0*
_output_shapes
: 

Mtraining/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdam/ReadVariableOp_2ReadVariableOpbeta_2*
dtype0*
_output_shapes
: 

Mtraining/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdam/ReadVariableOp_3ReadVariableOpepsilon*
dtype0*
_output_shapes
: 
¼
<training/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdamResourceApplyAdamin1_dense/kernel training/Adam/in1_dense/kernel/m training/Adam/in1_dense/kernel/v.training/Adam/Adam/update_in1_dense/kernel/Pow0training/Adam/Adam/update_in1_dense/kernel/Pow_1Ktraining/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdam/ReadVariableOpMtraining/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdam/ReadVariableOp_1Mtraining/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdam/ReadVariableOp_2Mtraining/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdam/ReadVariableOp_36training/Adam/gradients/in1_dense/MatMul_grad/MatMul_1*
T0*#
_class
loc:@in1_dense/kernel*
use_locking(

7training/Adam/Adam/update_in1_dense/bias/ReadVariableOpReadVariableOpiter"/device:CPU:0*
dtype0	*
_output_shapes
: 

.training/Adam/Adam/update_in1_dense/bias/add/yConst*
value	B	 R*!
_class
loc:@in1_dense/bias*
dtype0	*
_output_shapes
: 
ą
,training/Adam/Adam/update_in1_dense/bias/addAdd7training/Adam/Adam/update_in1_dense/bias/ReadVariableOp.training/Adam/Adam/update_in1_dense/bias/add/y*
T0	*!
_class
loc:@in1_dense/bias*
_output_shapes
: 
¶
-training/Adam/Adam/update_in1_dense/bias/CastCast,training/Adam/Adam/update_in1_dense/bias/add*

SrcT0	*!
_class
loc:@in1_dense/bias*
_output_shapes
: *

DstT0
z
;training/Adam/Adam/update_in1_dense/bias/Pow/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: 
ć
,training/Adam/Adam/update_in1_dense/bias/PowPow;training/Adam/Adam/update_in1_dense/bias/Pow/ReadVariableOp-training/Adam/Adam/update_in1_dense/bias/Cast*
T0*!
_class
loc:@in1_dense/bias*
_output_shapes
: 
|
=training/Adam/Adam/update_in1_dense/bias/Pow_1/ReadVariableOpReadVariableOpbeta_2*
dtype0*
_output_shapes
: 
ē
.training/Adam/Adam/update_in1_dense/bias/Pow_1Pow=training/Adam/Adam/update_in1_dense/bias/Pow_1/ReadVariableOp-training/Adam/Adam/update_in1_dense/bias/Cast*
_output_shapes
: *
T0*!
_class
loc:@in1_dense/bias

Itraining/Adam/Adam/update_in1_dense/bias/ResourceApplyAdam/ReadVariableOpReadVariableOplearning_rate*
dtype0*
_output_shapes
: 

Ktraining/Adam/Adam/update_in1_dense/bias/ResourceApplyAdam/ReadVariableOp_1ReadVariableOpbeta_1*
dtype0*
_output_shapes
: 

Ktraining/Adam/Adam/update_in1_dense/bias/ResourceApplyAdam/ReadVariableOp_2ReadVariableOpbeta_2*
dtype0*
_output_shapes
: 

Ktraining/Adam/Adam/update_in1_dense/bias/ResourceApplyAdam/ReadVariableOp_3ReadVariableOpepsilon*
dtype0*
_output_shapes
: 
Ŗ
:training/Adam/Adam/update_in1_dense/bias/ResourceApplyAdamResourceApplyAdamin1_dense/biastraining/Adam/in1_dense/bias/mtraining/Adam/in1_dense/bias/v,training/Adam/Adam/update_in1_dense/bias/Pow.training/Adam/Adam/update_in1_dense/bias/Pow_1Itraining/Adam/Adam/update_in1_dense/bias/ResourceApplyAdam/ReadVariableOpKtraining/Adam/Adam/update_in1_dense/bias/ResourceApplyAdam/ReadVariableOp_1Ktraining/Adam/Adam/update_in1_dense/bias/ResourceApplyAdam/ReadVariableOp_2Ktraining/Adam/Adam/update_in1_dense/bias/ResourceApplyAdam/ReadVariableOp_3:training/Adam/gradients/in1_dense/BiasAdd_grad/BiasAddGrad*
use_locking(*
T0*!
_class
loc:@in1_dense/bias

9training/Adam/Adam/update_in2_dense/kernel/ReadVariableOpReadVariableOpiter"/device:CPU:0*
dtype0	*
_output_shapes
: 

0training/Adam/Adam/update_in2_dense/kernel/add/yConst*
value	B	 R*#
_class
loc:@in2_dense/kernel*
dtype0	*
_output_shapes
: 
č
.training/Adam/Adam/update_in2_dense/kernel/addAdd9training/Adam/Adam/update_in2_dense/kernel/ReadVariableOp0training/Adam/Adam/update_in2_dense/kernel/add/y*
_output_shapes
: *
T0	*#
_class
loc:@in2_dense/kernel
¼
/training/Adam/Adam/update_in2_dense/kernel/CastCast.training/Adam/Adam/update_in2_dense/kernel/add*

SrcT0	*#
_class
loc:@in2_dense/kernel*
_output_shapes
: *

DstT0
|
=training/Adam/Adam/update_in2_dense/kernel/Pow/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: 
ė
.training/Adam/Adam/update_in2_dense/kernel/PowPow=training/Adam/Adam/update_in2_dense/kernel/Pow/ReadVariableOp/training/Adam/Adam/update_in2_dense/kernel/Cast*
T0*#
_class
loc:@in2_dense/kernel*
_output_shapes
: 
~
?training/Adam/Adam/update_in2_dense/kernel/Pow_1/ReadVariableOpReadVariableOpbeta_2*
dtype0*
_output_shapes
: 
ļ
0training/Adam/Adam/update_in2_dense/kernel/Pow_1Pow?training/Adam/Adam/update_in2_dense/kernel/Pow_1/ReadVariableOp/training/Adam/Adam/update_in2_dense/kernel/Cast*
T0*#
_class
loc:@in2_dense/kernel*
_output_shapes
: 

Ktraining/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdam/ReadVariableOpReadVariableOplearning_rate*
dtype0*
_output_shapes
: 

Mtraining/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdam/ReadVariableOp_1ReadVariableOpbeta_1*
dtype0*
_output_shapes
: 

Mtraining/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdam/ReadVariableOp_2ReadVariableOpbeta_2*
dtype0*
_output_shapes
: 

Mtraining/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdam/ReadVariableOp_3ReadVariableOpepsilon*
dtype0*
_output_shapes
: 
¼
<training/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdamResourceApplyAdamin2_dense/kernel training/Adam/in2_dense/kernel/m training/Adam/in2_dense/kernel/v.training/Adam/Adam/update_in2_dense/kernel/Pow0training/Adam/Adam/update_in2_dense/kernel/Pow_1Ktraining/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdam/ReadVariableOpMtraining/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdam/ReadVariableOp_1Mtraining/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdam/ReadVariableOp_2Mtraining/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdam/ReadVariableOp_36training/Adam/gradients/in2_dense/MatMul_grad/MatMul_1*
use_locking(*
T0*#
_class
loc:@in2_dense/kernel

7training/Adam/Adam/update_in2_dense/bias/ReadVariableOpReadVariableOpiter"/device:CPU:0*
dtype0	*
_output_shapes
: 

.training/Adam/Adam/update_in2_dense/bias/add/yConst*
value	B	 R*!
_class
loc:@in2_dense/bias*
dtype0	*
_output_shapes
: 
ą
,training/Adam/Adam/update_in2_dense/bias/addAdd7training/Adam/Adam/update_in2_dense/bias/ReadVariableOp.training/Adam/Adam/update_in2_dense/bias/add/y*
T0	*!
_class
loc:@in2_dense/bias*
_output_shapes
: 
¶
-training/Adam/Adam/update_in2_dense/bias/CastCast,training/Adam/Adam/update_in2_dense/bias/add*
_output_shapes
: *

DstT0*

SrcT0	*!
_class
loc:@in2_dense/bias
z
;training/Adam/Adam/update_in2_dense/bias/Pow/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: 
ć
,training/Adam/Adam/update_in2_dense/bias/PowPow;training/Adam/Adam/update_in2_dense/bias/Pow/ReadVariableOp-training/Adam/Adam/update_in2_dense/bias/Cast*
T0*!
_class
loc:@in2_dense/bias*
_output_shapes
: 
|
=training/Adam/Adam/update_in2_dense/bias/Pow_1/ReadVariableOpReadVariableOpbeta_2*
dtype0*
_output_shapes
: 
ē
.training/Adam/Adam/update_in2_dense/bias/Pow_1Pow=training/Adam/Adam/update_in2_dense/bias/Pow_1/ReadVariableOp-training/Adam/Adam/update_in2_dense/bias/Cast*
_output_shapes
: *
T0*!
_class
loc:@in2_dense/bias

Itraining/Adam/Adam/update_in2_dense/bias/ResourceApplyAdam/ReadVariableOpReadVariableOplearning_rate*
dtype0*
_output_shapes
: 

Ktraining/Adam/Adam/update_in2_dense/bias/ResourceApplyAdam/ReadVariableOp_1ReadVariableOpbeta_1*
dtype0*
_output_shapes
: 

Ktraining/Adam/Adam/update_in2_dense/bias/ResourceApplyAdam/ReadVariableOp_2ReadVariableOpbeta_2*
dtype0*
_output_shapes
: 

Ktraining/Adam/Adam/update_in2_dense/bias/ResourceApplyAdam/ReadVariableOp_3ReadVariableOpepsilon*
dtype0*
_output_shapes
: 
Ŗ
:training/Adam/Adam/update_in2_dense/bias/ResourceApplyAdamResourceApplyAdamin2_dense/biastraining/Adam/in2_dense/bias/mtraining/Adam/in2_dense/bias/v,training/Adam/Adam/update_in2_dense/bias/Pow.training/Adam/Adam/update_in2_dense/bias/Pow_1Itraining/Adam/Adam/update_in2_dense/bias/ResourceApplyAdam/ReadVariableOpKtraining/Adam/Adam/update_in2_dense/bias/ResourceApplyAdam/ReadVariableOp_1Ktraining/Adam/Adam/update_in2_dense/bias/ResourceApplyAdam/ReadVariableOp_2Ktraining/Adam/Adam/update_in2_dense/bias/ResourceApplyAdam/ReadVariableOp_3:training/Adam/gradients/in2_dense/BiasAdd_grad/BiasAddGrad*
use_locking(*
T0*!
_class
loc:@in2_dense/bias

6training/Adam/Adam/update_dense3/kernel/ReadVariableOpReadVariableOpiter"/device:CPU:0*
dtype0	*
_output_shapes
: 

-training/Adam/Adam/update_dense3/kernel/add/yConst*
value	B	 R* 
_class
loc:@dense3/kernel*
dtype0	*
_output_shapes
: 
Ü
+training/Adam/Adam/update_dense3/kernel/addAdd6training/Adam/Adam/update_dense3/kernel/ReadVariableOp-training/Adam/Adam/update_dense3/kernel/add/y*
T0	* 
_class
loc:@dense3/kernel*
_output_shapes
: 
³
,training/Adam/Adam/update_dense3/kernel/CastCast+training/Adam/Adam/update_dense3/kernel/add*

SrcT0	* 
_class
loc:@dense3/kernel*
_output_shapes
: *

DstT0
y
:training/Adam/Adam/update_dense3/kernel/Pow/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: 
ß
+training/Adam/Adam/update_dense3/kernel/PowPow:training/Adam/Adam/update_dense3/kernel/Pow/ReadVariableOp,training/Adam/Adam/update_dense3/kernel/Cast*
T0* 
_class
loc:@dense3/kernel*
_output_shapes
: 
{
<training/Adam/Adam/update_dense3/kernel/Pow_1/ReadVariableOpReadVariableOpbeta_2*
dtype0*
_output_shapes
: 
ć
-training/Adam/Adam/update_dense3/kernel/Pow_1Pow<training/Adam/Adam/update_dense3/kernel/Pow_1/ReadVariableOp,training/Adam/Adam/update_dense3/kernel/Cast*
T0* 
_class
loc:@dense3/kernel*
_output_shapes
: 

Htraining/Adam/Adam/update_dense3/kernel/ResourceApplyAdam/ReadVariableOpReadVariableOplearning_rate*
dtype0*
_output_shapes
: 

Jtraining/Adam/Adam/update_dense3/kernel/ResourceApplyAdam/ReadVariableOp_1ReadVariableOpbeta_1*
dtype0*
_output_shapes
: 

Jtraining/Adam/Adam/update_dense3/kernel/ResourceApplyAdam/ReadVariableOp_2ReadVariableOpbeta_2*
dtype0*
_output_shapes
: 

Jtraining/Adam/Adam/update_dense3/kernel/ResourceApplyAdam/ReadVariableOp_3ReadVariableOpepsilon*
dtype0*
_output_shapes
: 

9training/Adam/Adam/update_dense3/kernel/ResourceApplyAdamResourceApplyAdamdense3/kerneltraining/Adam/dense3/kernel/mtraining/Adam/dense3/kernel/v+training/Adam/Adam/update_dense3/kernel/Pow-training/Adam/Adam/update_dense3/kernel/Pow_1Htraining/Adam/Adam/update_dense3/kernel/ResourceApplyAdam/ReadVariableOpJtraining/Adam/Adam/update_dense3/kernel/ResourceApplyAdam/ReadVariableOp_1Jtraining/Adam/Adam/update_dense3/kernel/ResourceApplyAdam/ReadVariableOp_2Jtraining/Adam/Adam/update_dense3/kernel/ResourceApplyAdam/ReadVariableOp_33training/Adam/gradients/dense3/MatMul_grad/MatMul_1*
use_locking(*
T0* 
_class
loc:@dense3/kernel

4training/Adam/Adam/update_dense3/bias/ReadVariableOpReadVariableOpiter"/device:CPU:0*
dtype0	*
_output_shapes
: 

+training/Adam/Adam/update_dense3/bias/add/yConst*
value	B	 R*
_class
loc:@dense3/bias*
dtype0	*
_output_shapes
: 
Ō
)training/Adam/Adam/update_dense3/bias/addAdd4training/Adam/Adam/update_dense3/bias/ReadVariableOp+training/Adam/Adam/update_dense3/bias/add/y*
_output_shapes
: *
T0	*
_class
loc:@dense3/bias
­
*training/Adam/Adam/update_dense3/bias/CastCast)training/Adam/Adam/update_dense3/bias/add*
_output_shapes
: *

DstT0*

SrcT0	*
_class
loc:@dense3/bias
w
8training/Adam/Adam/update_dense3/bias/Pow/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: 
×
)training/Adam/Adam/update_dense3/bias/PowPow8training/Adam/Adam/update_dense3/bias/Pow/ReadVariableOp*training/Adam/Adam/update_dense3/bias/Cast*
T0*
_class
loc:@dense3/bias*
_output_shapes
: 
y
:training/Adam/Adam/update_dense3/bias/Pow_1/ReadVariableOpReadVariableOpbeta_2*
dtype0*
_output_shapes
: 
Ū
+training/Adam/Adam/update_dense3/bias/Pow_1Pow:training/Adam/Adam/update_dense3/bias/Pow_1/ReadVariableOp*training/Adam/Adam/update_dense3/bias/Cast*
T0*
_class
loc:@dense3/bias*
_output_shapes
: 

Ftraining/Adam/Adam/update_dense3/bias/ResourceApplyAdam/ReadVariableOpReadVariableOplearning_rate*
dtype0*
_output_shapes
: 

Htraining/Adam/Adam/update_dense3/bias/ResourceApplyAdam/ReadVariableOp_1ReadVariableOpbeta_1*
dtype0*
_output_shapes
: 

Htraining/Adam/Adam/update_dense3/bias/ResourceApplyAdam/ReadVariableOp_2ReadVariableOpbeta_2*
dtype0*
_output_shapes
: 

Htraining/Adam/Adam/update_dense3/bias/ResourceApplyAdam/ReadVariableOp_3ReadVariableOpepsilon*
dtype0*
_output_shapes
: 

7training/Adam/Adam/update_dense3/bias/ResourceApplyAdamResourceApplyAdamdense3/biastraining/Adam/dense3/bias/mtraining/Adam/dense3/bias/v)training/Adam/Adam/update_dense3/bias/Pow+training/Adam/Adam/update_dense3/bias/Pow_1Ftraining/Adam/Adam/update_dense3/bias/ResourceApplyAdam/ReadVariableOpHtraining/Adam/Adam/update_dense3/bias/ResourceApplyAdam/ReadVariableOp_1Htraining/Adam/Adam/update_dense3/bias/ResourceApplyAdam/ReadVariableOp_2Htraining/Adam/Adam/update_dense3/bias/ResourceApplyAdam/ReadVariableOp_37training/Adam/gradients/dense3/BiasAdd_grad/BiasAddGrad*
use_locking(*
T0*
_class
loc:@dense3/bias
Č
training/Adam/Adam/ConstConst8^training/Adam/Adam/update_dense3/bias/ResourceApplyAdam:^training/Adam/Adam/update_dense3/kernel/ResourceApplyAdam;^training/Adam/Adam/update_in1_dense/bias/ResourceApplyAdam=^training/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdam;^training/Adam/Adam/update_in2_dense/bias/ResourceApplyAdam=^training/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdam*
value	B	 R*
dtype0	*
_output_shapes
: 
j
&training/Adam/Adam/AssignAddVariableOpAssignAddVariableOpitertraining/Adam/Adam/Const*
dtype0	
õ
!training/Adam/Adam/ReadVariableOpReadVariableOpiter'^training/Adam/Adam/AssignAddVariableOp8^training/Adam/Adam/update_dense3/bias/ResourceApplyAdam:^training/Adam/Adam/update_dense3/kernel/ResourceApplyAdam;^training/Adam/Adam/update_in1_dense/bias/ResourceApplyAdam=^training/Adam/Adam/update_in1_dense/kernel/ResourceApplyAdam;^training/Adam/Adam/update_in2_dense/bias/ResourceApplyAdam=^training/Adam/Adam/update_in2_dense/kernel/ResourceApplyAdam*
dtype0	*
_output_shapes
: 
k
training_1/group_depsNoOp	^loss/mul^metrics/accuracy/Mean_2'^training/Adam/Adam/AssignAddVariableOp
W
Const_1Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_2Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_3Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_4Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_5Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_6Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_7Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
`
VarIsInitializedOpVarIsInitializedOptraining/Adam/in1_dense/bias/m*
_output_shapes
: 
R
VarIsInitializedOp_1VarIsInitializedOpin2_dense/bias*
_output_shapes
: 
_
VarIsInitializedOp_2VarIsInitializedOptraining/Adam/dense3/bias/m*
_output_shapes
: 
Q
VarIsInitializedOp_3VarIsInitializedOplearning_rate*
_output_shapes
: 
I
VarIsInitializedOp_4VarIsInitializedOptotal*
_output_shapes
: 
d
VarIsInitializedOp_5VarIsInitializedOp training/Adam/in2_dense/kernel/m*
_output_shapes
: 
d
VarIsInitializedOp_6VarIsInitializedOp training/Adam/in1_dense/kernel/m*
_output_shapes
: 
Q
VarIsInitializedOp_7VarIsInitializedOpdense3/kernel*
_output_shapes
: 
I
VarIsInitializedOp_8VarIsInitializedOpcount*
_output_shapes
: 
I
VarIsInitializedOp_9VarIsInitializedOpdecay*
_output_shapes
: 
b
VarIsInitializedOp_10VarIsInitializedOptraining/Adam/dense3/kernel/v*
_output_shapes
: 
c
VarIsInitializedOp_11VarIsInitializedOptraining/Adam/in2_dense/bias/v*
_output_shapes
: 
e
VarIsInitializedOp_12VarIsInitializedOp training/Adam/in2_dense/kernel/v*
_output_shapes
: 
U
VarIsInitializedOp_13VarIsInitializedOpin1_dense/kernel*
_output_shapes
: 
K
VarIsInitializedOp_14VarIsInitializedOpbeta_2*
_output_shapes
: 
`
VarIsInitializedOp_15VarIsInitializedOptraining/Adam/dense3/bias/v*
_output_shapes
: 
c
VarIsInitializedOp_16VarIsInitializedOptraining/Adam/in1_dense/bias/v*
_output_shapes
: 
S
VarIsInitializedOp_17VarIsInitializedOpin1_dense/bias*
_output_shapes
: 
L
VarIsInitializedOp_18VarIsInitializedOpepsilon*
_output_shapes
: 
e
VarIsInitializedOp_19VarIsInitializedOp training/Adam/in1_dense/kernel/v*
_output_shapes
: 
P
VarIsInitializedOp_20VarIsInitializedOpdense3/bias*
_output_shapes
: 
U
VarIsInitializedOp_21VarIsInitializedOpin2_dense/kernel*
_output_shapes
: 
I
VarIsInitializedOp_22VarIsInitializedOpiter*
_output_shapes
: 
b
VarIsInitializedOp_23VarIsInitializedOptraining/Adam/dense3/kernel/m*
_output_shapes
: 
K
VarIsInitializedOp_24VarIsInitializedOpbeta_1*
_output_shapes
: 
c
VarIsInitializedOp_25VarIsInitializedOptraining/Adam/in2_dense/bias/m*
_output_shapes
: 
ö
	init/NoOpNoOp^beta_1/Assign^beta_2/Assign^count/Assign^decay/Assign^dense3/bias/Assign^dense3/kernel/Assign^epsilon/Assign^in1_dense/bias/Assign^in1_dense/kernel/Assign^in2_dense/bias/Assign^in2_dense/kernel/Assign^learning_rate/Assign^total/Assign#^training/Adam/dense3/bias/m/Assign#^training/Adam/dense3/bias/v/Assign%^training/Adam/dense3/kernel/m/Assign%^training/Adam/dense3/kernel/v/Assign&^training/Adam/in1_dense/bias/m/Assign&^training/Adam/in1_dense/bias/v/Assign(^training/Adam/in1_dense/kernel/m/Assign(^training/Adam/in1_dense/kernel/v/Assign&^training/Adam/in2_dense/bias/m/Assign&^training/Adam/in2_dense/bias/v/Assign(^training/Adam/in2_dense/kernel/m/Assign(^training/Adam/in2_dense/kernel/v/Assign
0
init/NoOp_1NoOp^iter/Assign"/device:CPU:0
&
initNoOp
^init/NoOp^init/NoOp_1
W
Const_8Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
\
Const_9Const"/device:CPU:0*
valueB Bmodel*
dtype0*
_output_shapes
: 
X
Const_10Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_11Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_12Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_13Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
Ą
RestoreV2/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*g
value^B\BRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
r
RestoreV2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

	RestoreV2	RestoreV2Const_9RestoreV2/tensor_namesRestoreV2/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
B
IdentityIdentity	RestoreV2*
T0*
_output_shapes
:
]
AssignVariableOpAssignVariableOp training/Adam/in1_dense/kernel/mIdentity*
dtype0
Ā
RestoreV2_1/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*g
value^B\BRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
t
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_1	RestoreV2Const_9RestoreV2_1/tensor_namesRestoreV2_1/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_1IdentityRestoreV2_1*
_output_shapes
:*
T0
a
AssignVariableOp_1AssignVariableOp training/Adam/in1_dense/kernel/v
Identity_1*
dtype0
Ą
RestoreV2_2/tensor_namesConst"/device:CPU:0*e
value\BZBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_2	RestoreV2Const_9RestoreV2_2/tensor_namesRestoreV2_2/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_2IdentityRestoreV2_2*
T0*
_output_shapes
:
_
AssignVariableOp_2AssignVariableOptraining/Adam/in1_dense/bias/m
Identity_2*
dtype0
Ą
RestoreV2_3/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*e
value\BZBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
t
RestoreV2_3/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_3	RestoreV2Const_9RestoreV2_3/tensor_namesRestoreV2_3/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_3IdentityRestoreV2_3*
T0*
_output_shapes
:
_
AssignVariableOp_3AssignVariableOptraining/Adam/in1_dense/bias/v
Identity_3*
dtype0
X
Const_14Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
Ā
RestoreV2_4/tensor_namesConst"/device:CPU:0*g
value^B\BRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_4/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_4	RestoreV2Const_9RestoreV2_4/tensor_namesRestoreV2_4/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_4IdentityRestoreV2_4*
_output_shapes
:*
T0
a
AssignVariableOp_4AssignVariableOp training/Adam/in2_dense/kernel/m
Identity_4*
dtype0
Ā
RestoreV2_5/tensor_namesConst"/device:CPU:0*g
value^B\BRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_5/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_5	RestoreV2Const_9RestoreV2_5/tensor_namesRestoreV2_5/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_5IdentityRestoreV2_5*
T0*
_output_shapes
:
a
AssignVariableOp_5AssignVariableOp training/Adam/in2_dense/kernel/v
Identity_5*
dtype0
Ą
RestoreV2_6/tensor_namesConst"/device:CPU:0*e
value\BZBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_6/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_6	RestoreV2Const_9RestoreV2_6/tensor_namesRestoreV2_6/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_6IdentityRestoreV2_6*
T0*
_output_shapes
:
_
AssignVariableOp_6AssignVariableOptraining/Adam/in2_dense/bias/m
Identity_6*
dtype0
Ą
RestoreV2_7/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*e
value\BZBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
t
RestoreV2_7/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_7	RestoreV2Const_9RestoreV2_7/tensor_namesRestoreV2_7/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_7IdentityRestoreV2_7*
_output_shapes
:*
T0
_
AssignVariableOp_7AssignVariableOptraining/Adam/in2_dense/bias/v
Identity_7*
dtype0
X
Const_15Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_16Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
Ā
RestoreV2_8/tensor_namesConst"/device:CPU:0*g
value^B\BRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_8/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_8	RestoreV2Const_9RestoreV2_8/tensor_namesRestoreV2_8/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_8IdentityRestoreV2_8*
T0*
_output_shapes
:
^
AssignVariableOp_8AssignVariableOptraining/Adam/dense3/kernel/m
Identity_8*
dtype0
Ā
RestoreV2_9/tensor_namesConst"/device:CPU:0*g
value^B\BRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_9/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_9	RestoreV2Const_9RestoreV2_9/tensor_namesRestoreV2_9/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_9IdentityRestoreV2_9*
T0*
_output_shapes
:
^
AssignVariableOp_9AssignVariableOptraining/Adam/dense3/kernel/v
Identity_9*
dtype0
Į
RestoreV2_10/tensor_namesConst"/device:CPU:0*e
value\BZBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_10/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_10	RestoreV2Const_9RestoreV2_10/tensor_namesRestoreV2_10/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_10IdentityRestoreV2_10*
T0*
_output_shapes
:
^
AssignVariableOp_10AssignVariableOptraining/Adam/dense3/bias/mIdentity_10*
dtype0
Į
RestoreV2_11/tensor_namesConst"/device:CPU:0*e
value\BZBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_11/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_11	RestoreV2Const_9RestoreV2_11/tensor_namesRestoreV2_11/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_11IdentityRestoreV2_11*
_output_shapes
:*
T0
^
AssignVariableOp_11AssignVariableOptraining/Adam/dense3/bias/vIdentity_11*
dtype0
X
Const_17Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
§
RestoreV2_12/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_12/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_12	RestoreV2Const_9RestoreV2_12/tensor_namesRestoreV2_12/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_12IdentityRestoreV2_12*
T0*
_output_shapes
:
S
AssignVariableOp_12AssignVariableOpin1_dense/kernelIdentity_12*
dtype0
„
RestoreV2_13/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_13/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_13	RestoreV2Const_9RestoreV2_13/tensor_namesRestoreV2_13/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_13IdentityRestoreV2_13*
T0*
_output_shapes
:
Q
AssignVariableOp_13AssignVariableOpin1_dense/biasIdentity_13*
dtype0
§
RestoreV2_14/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_14/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_14	RestoreV2Const_9RestoreV2_14/tensor_namesRestoreV2_14/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
H
Identity_14IdentityRestoreV2_14*
_output_shapes
:*
T0
S
AssignVariableOp_14AssignVariableOpin2_dense/kernelIdentity_14*
dtype0
„
RestoreV2_15/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_15/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_15	RestoreV2Const_9RestoreV2_15/tensor_namesRestoreV2_15/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_15IdentityRestoreV2_15*
_output_shapes
:*
T0
Q
AssignVariableOp_15AssignVariableOpin2_dense/biasIdentity_15*
dtype0
§
RestoreV2_16/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_16/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_16	RestoreV2Const_9RestoreV2_16/tensor_namesRestoreV2_16/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_16IdentityRestoreV2_16*
T0*
_output_shapes
:
P
AssignVariableOp_16AssignVariableOpdense3/kernelIdentity_16*
dtype0
„
RestoreV2_17/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*I
value@B>B4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE
u
RestoreV2_17/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_17	RestoreV2Const_9RestoreV2_17/tensor_namesRestoreV2_17/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_17IdentityRestoreV2_17*
T0*
_output_shapes
:
N
AssignVariableOp_17AssignVariableOpdense3/biasIdentity_17*
dtype0

RestoreV2_18/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*>
value5B3B)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
u
RestoreV2_18/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_18	RestoreV2Const_9RestoreV2_18/tensor_namesRestoreV2_18/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2	
W
Identity_18IdentityRestoreV2_18"/device:CPU:0*
_output_shapes
:*
T0	
V
AssignVariableOp_18AssignVariableOpiterIdentity_18"/device:CPU:0*
dtype0	
£
RestoreV2_19/tensor_namesConst"/device:CPU:0*G
value>B<B2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_19/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_19	RestoreV2Const_9RestoreV2_19/tensor_namesRestoreV2_19/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_19IdentityRestoreV2_19*
T0*
_output_shapes
:
P
AssignVariableOp_19AssignVariableOplearning_rateIdentity_19*
dtype0

RestoreV2_20/tensor_namesConst"/device:CPU:0*?
value6B4B*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_20/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_20	RestoreV2Const_9RestoreV2_20/tensor_namesRestoreV2_20/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_20IdentityRestoreV2_20*
T0*
_output_shapes
:
H
AssignVariableOp_20AssignVariableOpdecayIdentity_20*
dtype0

RestoreV2_21/tensor_namesConst"/device:CPU:0*@
value7B5B+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_21/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_21	RestoreV2Const_9RestoreV2_21/tensor_namesRestoreV2_21/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_21IdentityRestoreV2_21*
T0*
_output_shapes
:
I
AssignVariableOp_21AssignVariableOpbeta_1Identity_21*
dtype0

RestoreV2_22/tensor_namesConst"/device:CPU:0*@
value7B5B+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_22/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_22	RestoreV2Const_9RestoreV2_22/tensor_namesRestoreV2_22/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_22IdentityRestoreV2_22*
T0*
_output_shapes
:
I
AssignVariableOp_22AssignVariableOpbeta_2Identity_22*
dtype0

RestoreV2_23/tensor_namesConst"/device:CPU:0*A
value8B6B,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_23/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_23	RestoreV2Const_9RestoreV2_23/tensor_namesRestoreV2_23/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
H
Identity_23IdentityRestoreV2_23*
T0*
_output_shapes
:
J
AssignVariableOp_23AssignVariableOpepsilonIdentity_23*
dtype0
X
Const_18Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
X
Const_19Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
­
SaveV2/tensor_namesConst"/device:CPU:0*Ö
valueĢBÉ!B/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-4/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-2/.ATTRIBUTES/OBJECT_CONFIG_JSONB(optimizer/.ATTRIBUTES/OBJECT_CONFIG_JSONB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH*
dtype0*
_output_shapes
:!
Æ
SaveV2/shape_and_slicesConst"/device:CPU:0*U
valueLBJ!B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:!
×	
SaveV2SaveV2Const_19SaveV2/tensor_namesSaveV2/shape_and_slicesConst_10Const_11Const_12Const_13Const_14Const_15Const_16Const_17$in1_dense/kernel/Read/ReadVariableOp"in1_dense/bias/Read/ReadVariableOp$in2_dense/kernel/Read/ReadVariableOp"in2_dense/bias/Read/ReadVariableOp!dense3/kernel/Read/ReadVariableOpdense3/bias/Read/ReadVariableOpiter/Read/ReadVariableOp!learning_rate/Read/ReadVariableOpdecay/Read/ReadVariableOpbeta_1/Read/ReadVariableOpbeta_2/Read/ReadVariableOpepsilon/Read/ReadVariableOp4training/Adam/in1_dense/kernel/m/Read/ReadVariableOp2training/Adam/in1_dense/bias/m/Read/ReadVariableOp4training/Adam/in2_dense/kernel/m/Read/ReadVariableOp2training/Adam/in2_dense/bias/m/Read/ReadVariableOp1training/Adam/dense3/kernel/m/Read/ReadVariableOp/training/Adam/dense3/bias/m/Read/ReadVariableOp4training/Adam/in1_dense/kernel/v/Read/ReadVariableOp2training/Adam/in1_dense/bias/v/Read/ReadVariableOp4training/Adam/in2_dense/kernel/v/Read/ReadVariableOp2training/Adam/in2_dense/bias/v/Read/ReadVariableOp1training/Adam/dense3/kernel/v/Read/ReadVariableOp/training/Adam/dense3/bias/v/Read/ReadVariableOpConst_18"/device:CPU:0*/
dtypes%
#2!	
Z
Identity_24IdentityConst_19^SaveV2"/device:CPU:0*
T0*
_output_shapes
: 
z
total_1/Initializer/zerosConst*
dtype0*
_output_shapes
: *
valueB
 *    *
_class
loc:@total_1
~
total_1VarHandleOp*
shape: *
shared_name	total_1*
_class
loc:@total_1*
dtype0*
_output_shapes
: 
_
(total_1/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal_1*
_output_shapes
: 
o
total_1/AssignAssignVariableOptotal_1total_1/Initializer/zeros*
_class
loc:@total_1*
dtype0
w
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_class
loc:@total_1*
dtype0*
_output_shapes
: 
z
count_1/Initializer/zerosConst*
valueB
 *    *
_class
loc:@count_1*
dtype0*
_output_shapes
: 
~
count_1VarHandleOp*
shared_name	count_1*
_class
loc:@count_1*
dtype0*
_output_shapes
: *
shape: 
_
(count_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount_1*
_output_shapes
: 
o
count_1/AssignAssignVariableOpcount_1count_1/Initializer/zeros*
_class
loc:@count_1*
dtype0
w
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_class
loc:@count_1*
dtype0*
_output_shapes
: 
K
Const_20Const*
valueB *
dtype0*
_output_shapes
: 
N
SumSummetrics/accuracy/Mean_2Const_20*
T0*
_output_shapes
: 
E
AssignAddVariableOpAssignAddVariableOptotal_1Sum*
dtype0
j
ReadVariableOpReadVariableOptotal_1^AssignAddVariableOp^Sum*
dtype0*
_output_shapes
: 
F
SizeConst*
value	B :*
dtype0*
_output_shapes
: 
B
CastCastSize*

SrcT0*
_output_shapes
: *

DstT0
^
AssignAddVariableOp_1AssignAddVariableOpcount_1Cast^AssignAddVariableOp*
dtype0
~
ReadVariableOp_1ReadVariableOpcount_1^AssignAddVariableOp^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
q
div_no_nan/ReadVariableOpReadVariableOptotal_1^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
s
div_no_nan/ReadVariableOp_1ReadVariableOpcount_1^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
o

div_no_nanDivNoNandiv_no_nan/ReadVariableOpdiv_no_nan/ReadVariableOp_1*
T0*
_output_shapes
: 
D
Identity_25Identity
div_no_nan*
_output_shapes
: *
T0
[
div_no_nan_1/ReadVariableOpReadVariableOptotal_1*
dtype0*
_output_shapes
: 
]
div_no_nan_1/ReadVariableOp_1ReadVariableOpcount_1*
dtype0*
_output_shapes
: 
u
div_no_nan_1DivNoNandiv_no_nan_1/ReadVariableOpdiv_no_nan_1/ReadVariableOp_1*
T0*
_output_shapes
: 
F
Identity_26Identitydiv_no_nan_1*
T0*
_output_shapes
: 
l
metric_op_wrapperConst^AssignAddVariableOp_1*
valueB *
dtype0*
_output_shapes
: 
Y
save/filename/inputConst*
valueB Bmodel*
dtype0*
_output_shapes
: 
n
save/filenamePlaceholderWithDefaultsave/filename/input*
dtype0*
_output_shapes
: *
shape: 
e

save/ConstPlaceholderWithDefaultsave/filename*
dtype0*
_output_shapes
: *
shape: 

save/Const_1Const*Ž
valueŌBŃ BŹ{"class_name": "Model", "config": {"input_layers": [["in1", 0, 0], ["in2", 0, 0]], "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 3], "dtype": "float32", "name": "in1", "sparse": false}, "inbound_nodes": [], "name": "in1"}, {"class_name": "InputLayer", "config": {"batch_input_shape": [null, 2], "dtype": "float32", "name": "in2", "sparse": false}, "inbound_nodes": [], "name": "in2"}, {"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in1_dense", "trainable": true, "units": 5, "use_bias": true}, "inbound_nodes": [["in1", 0, 0, {}]], "name": "in1_dense"}, {"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in2_dense", "trainable": true, "units": 10, "use_bias": true}, "inbound_nodes": [["in2", 0, 0, {}]], "name": "in2_dense"}, {"class_name": "Concatenate", "config": {"axis": -1, "dtype": "float32", "name": "merge", "trainable": true}, "inbound_nodes": [[["in1_dense", 0, 0, {}], ["in2_dense", 0, 0, {}]]], "name": "merge"}, {"class_name": "Dense", "config": {"activation": "sigmoid", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense3", "trainable": true, "units": 1, "use_bias": true}, "inbound_nodes": [["merge", 0, 0, {}]], "name": "dense3"}], "name": "model", "output_layers": ["dense3", 0, 0]}}*
dtype0*
_output_shapes
: 
Ģ
save/Const_2Const*
valueB B|{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 3], "dtype": "float32", "name": "in1", "sparse": false}}*
dtype0*
_output_shapes
: 
Ģ
save/Const_3Const*
dtype0*
_output_shapes
: *
valueB B|{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 2], "dtype": "float32", "name": "in2", "sparse": false}}
ŗ
save/Const_4Const*
dtype0*
_output_shapes
: *~
valueuBs Bm{"class_name": "Concatenate", "config": {"axis": -1, "dtype": "float32", "name": "merge", "trainable": true}}
ų
save/Const_5Const*»
value±B® B§{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in1_dense", "trainable": true, "units": 5, "use_bias": true}}*
dtype0*
_output_shapes
: 
ł
save/Const_6Const*¼
value²BÆ BØ{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in2_dense", "trainable": true, "units": 10, "use_bias": true}}*
dtype0*
_output_shapes
: 
ų
save/Const_7Const*»
value±B® B§{"class_name": "Dense", "config": {"activation": "sigmoid", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense3", "trainable": true, "units": 1, "use_bias": true}}*
dtype0*
_output_shapes
: 
N
save/VarIsInitializedOpVarIsInitializedOpcount_1*
_output_shapes
: 
P
save/VarIsInitializedOp_1VarIsInitializedOptotal_1*
_output_shapes
: 
3
	save/initNoOp^count_1/Assign^total_1/Assign
Ŗ
save/Const_8Const*
dtype0*
_output_shapes
: *ķ
valuećBą BŁ{"class_name": "Adam", "config": {"amsgrad": false, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "decay": 0.0, "epsilon": 1.0000000116860974e-07, "learning_rate": 0.0010000000474974513, "name": "Adam"}}

save/SaveV2/tensor_namesConst*
dtype0*
_output_shapes
: *ø
value®B« B/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-4/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-2/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB(optimizer/.ATTRIBUTES/OBJECT_CONFIG_JSONB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE
£
save/SaveV2/shape_and_slicesConst*
dtype0*
_output_shapes
: *S
valueJBH B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 
ī	
save/SaveV2SaveV2
save/Constsave/SaveV2/tensor_namessave/SaveV2/shape_and_slicessave/Const_1save/Const_2save/Const_3save/Const_4save/Const_5"in1_dense/bias/Read/ReadVariableOp2training/Adam/in1_dense/bias/m/Read/ReadVariableOp2training/Adam/in1_dense/bias/v/Read/ReadVariableOp$in1_dense/kernel/Read/ReadVariableOp4training/Adam/in1_dense/kernel/m/Read/ReadVariableOp4training/Adam/in1_dense/kernel/v/Read/ReadVariableOpsave/Const_6"in2_dense/bias/Read/ReadVariableOp2training/Adam/in2_dense/bias/m/Read/ReadVariableOp2training/Adam/in2_dense/bias/v/Read/ReadVariableOp$in2_dense/kernel/Read/ReadVariableOp4training/Adam/in2_dense/kernel/m/Read/ReadVariableOp4training/Adam/in2_dense/kernel/v/Read/ReadVariableOpsave/Const_7dense3/bias/Read/ReadVariableOp/training/Adam/dense3/bias/m/Read/ReadVariableOp/training/Adam/dense3/bias/v/Read/ReadVariableOp!dense3/kernel/Read/ReadVariableOp1training/Adam/dense3/kernel/m/Read/ReadVariableOp1training/Adam/dense3/kernel/v/Read/ReadVariableOpsave/Const_8beta_1/Read/ReadVariableOpbeta_2/Read/ReadVariableOpdecay/Read/ReadVariableOpepsilon/Read/ReadVariableOpiter/Read/ReadVariableOp!learning_rate/Read/ReadVariableOp*.
dtypes$
"2 	
}
save/control_dependencyIdentity
save/Const^save/SaveV2*
T0*
_class
loc:@save/Const*
_output_shapes
: 

save/RestoreV2/tensor_namesConst"/device:CPU:0*ø
value®B« B/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-4/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-2/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-2/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-2/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB(optimizer/.ATTRIBUTES/OBJECT_CONFIG_JSONB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
: 
µ
save/RestoreV2/shape_and_slicesConst"/device:CPU:0*S
valueJBH B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
: 
½
save/RestoreV2	RestoreV2
save/Constsave/RestoreV2/tensor_namessave/RestoreV2/shape_and_slices"/device:CPU:0*
_output_shapes
::::::::::::::::::::::::::::::::*.
dtypes$
"2 	

	save/NoOpNoOp

save/NoOp_1NoOp

save/NoOp_2NoOp

save/NoOp_3NoOp

save/NoOp_4NoOp
N
save/IdentityIdentitysave/RestoreV2:5*
T0*
_output_shapes
:
U
save/AssignVariableOpAssignVariableOpin1_dense/biassave/Identity*
dtype0
P
save/Identity_1Identitysave/RestoreV2:6*
T0*
_output_shapes
:
i
save/AssignVariableOp_1AssignVariableOptraining/Adam/in1_dense/bias/msave/Identity_1*
dtype0
P
save/Identity_2Identitysave/RestoreV2:7*
T0*
_output_shapes
:
i
save/AssignVariableOp_2AssignVariableOptraining/Adam/in1_dense/bias/vsave/Identity_2*
dtype0
P
save/Identity_3Identitysave/RestoreV2:8*
_output_shapes
:*
T0
[
save/AssignVariableOp_3AssignVariableOpin1_dense/kernelsave/Identity_3*
dtype0
P
save/Identity_4Identitysave/RestoreV2:9*
_output_shapes
:*
T0
k
save/AssignVariableOp_4AssignVariableOp training/Adam/in1_dense/kernel/msave/Identity_4*
dtype0
Q
save/Identity_5Identitysave/RestoreV2:10*
T0*
_output_shapes
:
k
save/AssignVariableOp_5AssignVariableOp training/Adam/in1_dense/kernel/vsave/Identity_5*
dtype0

save/NoOp_5NoOp
Q
save/Identity_6Identitysave/RestoreV2:12*
T0*
_output_shapes
:
Y
save/AssignVariableOp_6AssignVariableOpin2_dense/biassave/Identity_6*
dtype0
Q
save/Identity_7Identitysave/RestoreV2:13*
_output_shapes
:*
T0
i
save/AssignVariableOp_7AssignVariableOptraining/Adam/in2_dense/bias/msave/Identity_7*
dtype0
Q
save/Identity_8Identitysave/RestoreV2:14*
T0*
_output_shapes
:
i
save/AssignVariableOp_8AssignVariableOptraining/Adam/in2_dense/bias/vsave/Identity_8*
dtype0
Q
save/Identity_9Identitysave/RestoreV2:15*
_output_shapes
:*
T0
[
save/AssignVariableOp_9AssignVariableOpin2_dense/kernelsave/Identity_9*
dtype0
R
save/Identity_10Identitysave/RestoreV2:16*
T0*
_output_shapes
:
m
save/AssignVariableOp_10AssignVariableOp training/Adam/in2_dense/kernel/msave/Identity_10*
dtype0
R
save/Identity_11Identitysave/RestoreV2:17*
T0*
_output_shapes
:
m
save/AssignVariableOp_11AssignVariableOp training/Adam/in2_dense/kernel/vsave/Identity_11*
dtype0

save/NoOp_6NoOp
R
save/Identity_12Identitysave/RestoreV2:19*
T0*
_output_shapes
:
X
save/AssignVariableOp_12AssignVariableOpdense3/biassave/Identity_12*
dtype0
R
save/Identity_13Identitysave/RestoreV2:20*
T0*
_output_shapes
:
h
save/AssignVariableOp_13AssignVariableOptraining/Adam/dense3/bias/msave/Identity_13*
dtype0
R
save/Identity_14Identitysave/RestoreV2:21*
T0*
_output_shapes
:
h
save/AssignVariableOp_14AssignVariableOptraining/Adam/dense3/bias/vsave/Identity_14*
dtype0
R
save/Identity_15Identitysave/RestoreV2:22*
_output_shapes
:*
T0
Z
save/AssignVariableOp_15AssignVariableOpdense3/kernelsave/Identity_15*
dtype0
R
save/Identity_16Identitysave/RestoreV2:23*
T0*
_output_shapes
:
j
save/AssignVariableOp_16AssignVariableOptraining/Adam/dense3/kernel/msave/Identity_16*
dtype0
R
save/Identity_17Identitysave/RestoreV2:24*
_output_shapes
:*
T0
j
save/AssignVariableOp_17AssignVariableOptraining/Adam/dense3/kernel/vsave/Identity_17*
dtype0

save/NoOp_7NoOp
R
save/Identity_18Identitysave/RestoreV2:26*
T0*
_output_shapes
:
S
save/AssignVariableOp_18AssignVariableOpbeta_1save/Identity_18*
dtype0
R
save/Identity_19Identitysave/RestoreV2:27*
T0*
_output_shapes
:
S
save/AssignVariableOp_19AssignVariableOpbeta_2save/Identity_19*
dtype0
R
save/Identity_20Identitysave/RestoreV2:28*
T0*
_output_shapes
:
R
save/AssignVariableOp_20AssignVariableOpdecaysave/Identity_20*
dtype0
R
save/Identity_21Identitysave/RestoreV2:29*
_output_shapes
:*
T0
T
save/AssignVariableOp_21AssignVariableOpepsilonsave/Identity_21*
dtype0
a
save/Identity_22Identitysave/RestoreV2:30"/device:CPU:0*
T0	*
_output_shapes
:
`
save/AssignVariableOp_22AssignVariableOpitersave/Identity_22"/device:CPU:0*
dtype0	
R
save/Identity_23Identitysave/RestoreV2:31*
_output_shapes
:*
T0
Z
save/AssignVariableOp_23AssignVariableOplearning_ratesave/Identity_23*
dtype0
ģ
save/restore_all/NoOpNoOp^save/AssignVariableOp^save/AssignVariableOp_1^save/AssignVariableOp_10^save/AssignVariableOp_11^save/AssignVariableOp_12^save/AssignVariableOp_13^save/AssignVariableOp_14^save/AssignVariableOp_15^save/AssignVariableOp_16^save/AssignVariableOp_17^save/AssignVariableOp_18^save/AssignVariableOp_19^save/AssignVariableOp_2^save/AssignVariableOp_20^save/AssignVariableOp_21^save/AssignVariableOp_23^save/AssignVariableOp_3^save/AssignVariableOp_4^save/AssignVariableOp_5^save/AssignVariableOp_6^save/AssignVariableOp_7^save/AssignVariableOp_8^save/AssignVariableOp_9
^save/NoOp^save/NoOp_1^save/NoOp_2^save/NoOp_3^save/NoOp_4^save/NoOp_5^save/NoOp_6^save/NoOp_7
I
save/restore_all/NoOp_1NoOp^save/AssignVariableOp_22"/device:CPU:0
J
save/restore_allNoOp^save/restore_all/NoOp^save/restore_all/NoOp_1
0
init_1NoOp^count_1/Assign^total_1/Assign"D
save/Const:0save/control_dependency:0save/restore_all 5 @F8"
trainable_variablesöó

in1_dense/kernel:0in1_dense/kernel/Assign&in1_dense/kernel/Read/ReadVariableOp:0(2-in1_dense/kernel/Initializer/random_uniform:08
w
in1_dense/bias:0in1_dense/bias/Assign$in1_dense/bias/Read/ReadVariableOp:0(2"in1_dense/bias/Initializer/zeros:08

in2_dense/kernel:0in2_dense/kernel/Assign&in2_dense/kernel/Read/ReadVariableOp:0(2-in2_dense/kernel/Initializer/random_uniform:08
w
in2_dense/bias:0in2_dense/bias/Assign$in2_dense/bias/Read/ReadVariableOp:0(2"in2_dense/bias/Initializer/zeros:08
|
dense3/kernel:0dense3/kernel/Assign#dense3/kernel/Read/ReadVariableOp:0(2*dense3/kernel/Initializer/random_uniform:08
k
dense3/bias:0dense3/bias/Assign!dense3/bias/Read/ReadVariableOp:0(2dense3/bias/Initializer/zeros:08"Ķ
local_variables¹¶
Y
	total_1:0total_1/Assigntotal_1/Read/ReadVariableOp:0(2total_1/Initializer/zeros:0
Y
	count_1:0count_1/Assigncount_1/Read/ReadVariableOp:0(2count_1/Initializer/zeros:0"ö
cond_contextäą
Ś

qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/cond_textqloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_t:0 *ū
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:0
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1
qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_t:0ę
qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0Ś
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1
īq
sloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/cond_text_1qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_f:0*Ė5
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Merge:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Merge:1
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:1
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:1
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:0
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:2
¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1
£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims:0
Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
Ŗloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1
„loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim:0
”loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1:0
 loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat:0
„loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims:0
¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const:0
¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Shape:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0
qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_f:0
floss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank:0
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1:0
floss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch:0
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0ę
qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:02ł+
ö+
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/cond_textloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0 *Ę(
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:0
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:2
¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1
£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims:0
Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
Ŗloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1
„loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim:0
”loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1:0
 loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat:0
„loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims:0
¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const:0
¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Shape:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0Ņ
¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0Ŗloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1Ö
Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:02Å
Ā
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/cond_text_1loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0*

loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:1
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0¢
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0

nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/cond_textnloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t:0 *æ
yloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency:0
nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0
oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t:0ą
nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0
÷
ploss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/cond_text_1nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f:0*”
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch:0
vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1:0
vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2:0
vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7:0
{loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1:0
nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0
oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f:0
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0
oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0ē
oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge:0tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch:0ą
nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0ā
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1:0Ž
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3:0į
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2:0"`
global_stepQO
M
iter:0iter/Assigniter/Read/ReadVariableOp:0(2iter/Initializer/zeros:0"ź
	variablesÜŁ

in1_dense/kernel:0in1_dense/kernel/Assign&in1_dense/kernel/Read/ReadVariableOp:0(2-in1_dense/kernel/Initializer/random_uniform:08
w
in1_dense/bias:0in1_dense/bias/Assign$in1_dense/bias/Read/ReadVariableOp:0(2"in1_dense/bias/Initializer/zeros:08

in2_dense/kernel:0in2_dense/kernel/Assign&in2_dense/kernel/Read/ReadVariableOp:0(2-in2_dense/kernel/Initializer/random_uniform:08
w
in2_dense/bias:0in2_dense/bias/Assign$in2_dense/bias/Read/ReadVariableOp:0(2"in2_dense/bias/Initializer/zeros:08
|
dense3/kernel:0dense3/kernel/Assign#dense3/kernel/Read/ReadVariableOp:0(2*dense3/kernel/Initializer/random_uniform:08
k
dense3/bias:0dense3/bias/Assign!dense3/bias/Read/ReadVariableOp:0(2dense3/bias/Initializer/zeros:08
M
iter:0iter/Assigniter/Read/ReadVariableOp:0(2iter/Initializer/zeros:0
y
learning_rate:0learning_rate/Assign#learning_rate/Read/ReadVariableOp:0(2)learning_rate/Initializer/initial_value:0
Y
decay:0decay/Assigndecay/Read/ReadVariableOp:0(2!decay/Initializer/initial_value:0
]
beta_1:0beta_1/Assignbeta_1/Read/ReadVariableOp:0(2"beta_1/Initializer/initial_value:0
]
beta_2:0beta_2/Assignbeta_2/Read/ReadVariableOp:0(2"beta_2/Initializer/initial_value:0
a
	epsilon:0epsilon/Assignepsilon/Read/ReadVariableOp:0(2#epsilon/Initializer/initial_value:0
½
"training/Adam/in1_dense/kernel/m:0'training/Adam/in1_dense/kernel/m/Assign6training/Adam/in1_dense/kernel/m/Read/ReadVariableOp:0(24training/Adam/in1_dense/kernel/m/Initializer/zeros:0
µ
 training/Adam/in1_dense/bias/m:0%training/Adam/in1_dense/bias/m/Assign4training/Adam/in1_dense/bias/m/Read/ReadVariableOp:0(22training/Adam/in1_dense/bias/m/Initializer/zeros:0
½
"training/Adam/in2_dense/kernel/m:0'training/Adam/in2_dense/kernel/m/Assign6training/Adam/in2_dense/kernel/m/Read/ReadVariableOp:0(24training/Adam/in2_dense/kernel/m/Initializer/zeros:0
µ
 training/Adam/in2_dense/bias/m:0%training/Adam/in2_dense/bias/m/Assign4training/Adam/in2_dense/bias/m/Read/ReadVariableOp:0(22training/Adam/in2_dense/bias/m/Initializer/zeros:0
±
training/Adam/dense3/kernel/m:0$training/Adam/dense3/kernel/m/Assign3training/Adam/dense3/kernel/m/Read/ReadVariableOp:0(21training/Adam/dense3/kernel/m/Initializer/zeros:0
©
training/Adam/dense3/bias/m:0"training/Adam/dense3/bias/m/Assign1training/Adam/dense3/bias/m/Read/ReadVariableOp:0(2/training/Adam/dense3/bias/m/Initializer/zeros:0
½
"training/Adam/in1_dense/kernel/v:0'training/Adam/in1_dense/kernel/v/Assign6training/Adam/in1_dense/kernel/v/Read/ReadVariableOp:0(24training/Adam/in1_dense/kernel/v/Initializer/zeros:0
µ
 training/Adam/in1_dense/bias/v:0%training/Adam/in1_dense/bias/v/Assign4training/Adam/in1_dense/bias/v/Read/ReadVariableOp:0(22training/Adam/in1_dense/bias/v/Initializer/zeros:0
½
"training/Adam/in2_dense/kernel/v:0'training/Adam/in2_dense/kernel/v/Assign6training/Adam/in2_dense/kernel/v/Read/ReadVariableOp:0(24training/Adam/in2_dense/kernel/v/Initializer/zeros:0
µ
 training/Adam/in2_dense/bias/v:0%training/Adam/in2_dense/bias/v/Assign4training/Adam/in2_dense/bias/v/Read/ReadVariableOp:0(22training/Adam/in2_dense/bias/v/Initializer/zeros:0
±
training/Adam/dense3/kernel/v:0$training/Adam/dense3/kernel/v/Assign3training/Adam/dense3/kernel/v/Read/ReadVariableOp:0(21training/Adam/dense3/kernel/v/Initializer/zeros:0
©
training/Adam/dense3/bias/v:0"training/Adam/dense3/bias/v/Assign1training/Adam/dense3/bias/v/Read/ReadVariableOp:0(2/training/Adam/dense3/bias/v/Initializer/zeros:0*Q
__saved_model_train_op75
__saved_model_train_op
training_1/group_deps*ł
trainļ
@
dense3_target/
dense3_target:0’’’’’’’’’’’’’’’’’’
#
in1
in1:0’’’’’’’’’
#
in2
in2:0’’’’’’’’’9
metrics/accuracy/update_op
metric_op_wrapper:0 -
metrics/accuracy/value
Identity_26:0 =
predictions/dense3'
dense3/Sigmoid:0’’’’’’’’’
loss

loss/mul:0 tensorflow/supervised/training*@
__saved_model_init_op'%
__saved_model_init_op
init_1Éś
¦ó
:
Add
x"T
y"T
z"T"
Ttype:
2	
P
Assert
	condition
	
data2T"
T
list(type)(0"
	summarizeint
E
AssignAddVariableOp
resource
value"dtype"
dtypetype
B
AssignVariableOp
resource
value"dtype"
dtypetype
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
N
Cast	
x"SrcT	
y"DstT"
SrcTtype"
DstTtype"
Truncatebool( 
h
ConcatV2
values"T*N
axis"Tidx
output"T"
Nint(0"	
Ttype"
Tidxtype0:
2	
8
Const
output"dtype"
valuetensor"
dtypetype
¹
DenseToDenseSetOperation	
set1"T	
set2"T
result_indices	
result_values"T
result_shape	"
set_operationstring"
validate_indicesbool("
Ttype:
	2	
5
DivNoNan
x"T
y"T
z"T"
Ttype:
2
B
Equal
x"T
y"T
z
"
Ttype:
2	

W

ExpandDims

input"T
dim"Tdim
output"T"	
Ttype"
Tdimtype0:
2	
^
Fill
dims"
index_type

value"T
output"T"	
Ttype"

index_typetype0:
2	
=
Greater
x"T
y"T
z
"
Ttype:
2	
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	

Mean

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
N
Merge
inputs"T*N
output"T
value_index"	
Ttype"
Nint(0
=
Mul
x"T
y"T
z"T"
Ttype:
2	

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
X
PlaceholderWithDefault
input"dtype
output"dtype"
dtypetype"
shapeshape
~
RandomUniform

shape"T
output"dtype"
seedint "
seed2int "
dtypetype:
2"
Ttype:
2	
@
ReadVariableOp
resource
value"dtype"
dtypetype
E
Relu
features"T
activations"T"
Ttype:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
P
Shape

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
0
Sigmoid
x"T
y"T"
Ttype:

2
O
Size

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
G
SquaredDifference
x"T
y"T
z"T"
Ttype:

2	
:
Sub
x"T
y"T
z"T"
Ttype:
2	

Sum

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
M
Switch	
data"T
pred

output_false"T
output_true"T"	
Ttype
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape
9
VarIsInitializedOp
resource
is_initialized
"eval*2.0.0-alpha02v1.12.0-9492-g2c319fb4158¬
f
in1Placeholder*
dtype0*'
_output_shapes
:’’’’’’’’’*
shape:’’’’’’’’’
f
in2Placeholder*
dtype0*'
_output_shapes
:’’’’’’’’’*
shape:’’’’’’’’’
§
1in1_dense/kernel/Initializer/random_uniform/shapeConst*
valueB"      *#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes
:

/in1_dense/kernel/Initializer/random_uniform/minConst*
valueB
 *×³]æ*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes
: 

/in1_dense/kernel/Initializer/random_uniform/maxConst*
valueB
 *×³]?*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes
: 
Ų
9in1_dense/kernel/Initializer/random_uniform/RandomUniformRandomUniform1in1_dense/kernel/Initializer/random_uniform/shape*
T0*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes

:
Ž
/in1_dense/kernel/Initializer/random_uniform/subSub/in1_dense/kernel/Initializer/random_uniform/max/in1_dense/kernel/Initializer/random_uniform/min*
T0*#
_class
loc:@in1_dense/kernel*
_output_shapes
: 
š
/in1_dense/kernel/Initializer/random_uniform/mulMul9in1_dense/kernel/Initializer/random_uniform/RandomUniform/in1_dense/kernel/Initializer/random_uniform/sub*
_output_shapes

:*
T0*#
_class
loc:@in1_dense/kernel
ā
+in1_dense/kernel/Initializer/random_uniformAdd/in1_dense/kernel/Initializer/random_uniform/mul/in1_dense/kernel/Initializer/random_uniform/min*
T0*#
_class
loc:@in1_dense/kernel*
_output_shapes

:
”
in1_dense/kernelVarHandleOp*
shape
:*!
shared_namein1_dense/kernel*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes
: 
q
1in1_dense/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpin1_dense/kernel*
_output_shapes
: 

in1_dense/kernel/AssignAssignVariableOpin1_dense/kernel+in1_dense/kernel/Initializer/random_uniform*#
_class
loc:@in1_dense/kernel*
dtype0

$in1_dense/kernel/Read/ReadVariableOpReadVariableOpin1_dense/kernel*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes

:

 in1_dense/bias/Initializer/zerosConst*
valueB*    *!
_class
loc:@in1_dense/bias*
dtype0*
_output_shapes
:

in1_dense/biasVarHandleOp*!
_class
loc:@in1_dense/bias*
dtype0*
_output_shapes
: *
shape:*
shared_namein1_dense/bias
m
/in1_dense/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpin1_dense/bias*
_output_shapes
: 

in1_dense/bias/AssignAssignVariableOpin1_dense/bias in1_dense/bias/Initializer/zeros*!
_class
loc:@in1_dense/bias*
dtype0

"in1_dense/bias/Read/ReadVariableOpReadVariableOpin1_dense/bias*!
_class
loc:@in1_dense/bias*
dtype0*
_output_shapes
:
p
in1_dense/MatMul/ReadVariableOpReadVariableOpin1_dense/kernel*
dtype0*
_output_shapes

:
r
in1_dense/MatMulMatMulin1in1_dense/MatMul/ReadVariableOp*'
_output_shapes
:’’’’’’’’’*
T0
k
 in1_dense/BiasAdd/ReadVariableOpReadVariableOpin1_dense/bias*
dtype0*
_output_shapes
:

in1_dense/BiasAddBiasAddin1_dense/MatMul in1_dense/BiasAdd/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’
[
in1_dense/ReluReluin1_dense/BiasAdd*
T0*'
_output_shapes
:’’’’’’’’’
§
1in2_dense/kernel/Initializer/random_uniform/shapeConst*
dtype0*
_output_shapes
:*
valueB"   
   *#
_class
loc:@in2_dense/kernel

/in2_dense/kernel/Initializer/random_uniform/minConst*
valueB
 *ó5æ*#
_class
loc:@in2_dense/kernel*
dtype0*
_output_shapes
: 

/in2_dense/kernel/Initializer/random_uniform/maxConst*
dtype0*
_output_shapes
: *
valueB
 *ó5?*#
_class
loc:@in2_dense/kernel
Ų
9in2_dense/kernel/Initializer/random_uniform/RandomUniformRandomUniform1in2_dense/kernel/Initializer/random_uniform/shape*
T0*#
_class
loc:@in2_dense/kernel*
dtype0*
_output_shapes

:

Ž
/in2_dense/kernel/Initializer/random_uniform/subSub/in2_dense/kernel/Initializer/random_uniform/max/in2_dense/kernel/Initializer/random_uniform/min*
T0*#
_class
loc:@in2_dense/kernel*
_output_shapes
: 
š
/in2_dense/kernel/Initializer/random_uniform/mulMul9in2_dense/kernel/Initializer/random_uniform/RandomUniform/in2_dense/kernel/Initializer/random_uniform/sub*
T0*#
_class
loc:@in2_dense/kernel*
_output_shapes

:

ā
+in2_dense/kernel/Initializer/random_uniformAdd/in2_dense/kernel/Initializer/random_uniform/mul/in2_dense/kernel/Initializer/random_uniform/min*
T0*#
_class
loc:@in2_dense/kernel*
_output_shapes

:

”
in2_dense/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape
:
*!
shared_namein2_dense/kernel*#
_class
loc:@in2_dense/kernel
q
1in2_dense/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpin2_dense/kernel*
_output_shapes
: 

in2_dense/kernel/AssignAssignVariableOpin2_dense/kernel+in2_dense/kernel/Initializer/random_uniform*
dtype0*#
_class
loc:@in2_dense/kernel

$in2_dense/kernel/Read/ReadVariableOpReadVariableOpin2_dense/kernel*#
_class
loc:@in2_dense/kernel*
dtype0*
_output_shapes

:


 in2_dense/bias/Initializer/zerosConst*
valueB
*    *!
_class
loc:@in2_dense/bias*
dtype0*
_output_shapes
:


in2_dense/biasVarHandleOp*
shared_namein2_dense/bias*!
_class
loc:@in2_dense/bias*
dtype0*
_output_shapes
: *
shape:

m
/in2_dense/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpin2_dense/bias*
_output_shapes
: 

in2_dense/bias/AssignAssignVariableOpin2_dense/bias in2_dense/bias/Initializer/zeros*!
_class
loc:@in2_dense/bias*
dtype0

"in2_dense/bias/Read/ReadVariableOpReadVariableOpin2_dense/bias*!
_class
loc:@in2_dense/bias*
dtype0*
_output_shapes
:

p
in2_dense/MatMul/ReadVariableOpReadVariableOpin2_dense/kernel*
dtype0*
_output_shapes

:

r
in2_dense/MatMulMatMulin2in2_dense/MatMul/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’

k
 in2_dense/BiasAdd/ReadVariableOpReadVariableOpin2_dense/bias*
dtype0*
_output_shapes
:


in2_dense/BiasAddBiasAddin2_dense/MatMul in2_dense/BiasAdd/ReadVariableOp*'
_output_shapes
:’’’’’’’’’
*
T0
[
in2_dense/ReluReluin2_dense/BiasAdd*
T0*'
_output_shapes
:’’’’’’’’’

S
merge/concat/axisConst*
dtype0*
_output_shapes
: *
value	B :

merge/concatConcatV2in1_dense/Reluin2_dense/Relumerge/concat/axis*
T0*
N*'
_output_shapes
:’’’’’’’’’
”
.dense3/kernel/Initializer/random_uniform/shapeConst*
valueB"      * 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes
:

,dense3/kernel/Initializer/random_uniform/minConst*
valueB
 *qÄæ* 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes
: 

,dense3/kernel/Initializer/random_uniform/maxConst*
dtype0*
_output_shapes
: *
valueB
 *qÄ?* 
_class
loc:@dense3/kernel
Ļ
6dense3/kernel/Initializer/random_uniform/RandomUniformRandomUniform.dense3/kernel/Initializer/random_uniform/shape*
T0* 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes

:
Ņ
,dense3/kernel/Initializer/random_uniform/subSub,dense3/kernel/Initializer/random_uniform/max,dense3/kernel/Initializer/random_uniform/min*
_output_shapes
: *
T0* 
_class
loc:@dense3/kernel
ä
,dense3/kernel/Initializer/random_uniform/mulMul6dense3/kernel/Initializer/random_uniform/RandomUniform,dense3/kernel/Initializer/random_uniform/sub*
T0* 
_class
loc:@dense3/kernel*
_output_shapes

:
Ö
(dense3/kernel/Initializer/random_uniformAdd,dense3/kernel/Initializer/random_uniform/mul,dense3/kernel/Initializer/random_uniform/min*
T0* 
_class
loc:@dense3/kernel*
_output_shapes

:

dense3/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape
:*
shared_namedense3/kernel* 
_class
loc:@dense3/kernel
k
.dense3/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense3/kernel*
_output_shapes
: 

dense3/kernel/AssignAssignVariableOpdense3/kernel(dense3/kernel/Initializer/random_uniform* 
_class
loc:@dense3/kernel*
dtype0

!dense3/kernel/Read/ReadVariableOpReadVariableOpdense3/kernel*
dtype0*
_output_shapes

:* 
_class
loc:@dense3/kernel

dense3/bias/Initializer/zerosConst*
valueB*    *
_class
loc:@dense3/bias*
dtype0*
_output_shapes
:

dense3/biasVarHandleOp*
dtype0*
_output_shapes
: *
shape:*
shared_namedense3/bias*
_class
loc:@dense3/bias
g
,dense3/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense3/bias*
_output_shapes
: 

dense3/bias/AssignAssignVariableOpdense3/biasdense3/bias/Initializer/zeros*
dtype0*
_class
loc:@dense3/bias

dense3/bias/Read/ReadVariableOpReadVariableOpdense3/bias*
_class
loc:@dense3/bias*
dtype0*
_output_shapes
:
j
dense3/MatMul/ReadVariableOpReadVariableOpdense3/kernel*
dtype0*
_output_shapes

:
u
dense3/MatMulMatMulmerge/concatdense3/MatMul/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’
e
dense3/BiasAdd/ReadVariableOpReadVariableOpdense3/bias*
dtype0*
_output_shapes
:
y
dense3/BiasAddBiasAdddense3/MatMuldense3/BiasAdd/ReadVariableOp*'
_output_shapes
:’’’’’’’’’*
T0
[
dense3/SigmoidSigmoiddense3/BiasAdd*'
_output_shapes
:’’’’’’’’’*
T0

dense3_targetPlaceholder*
dtype0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’*%
shape:’’’’’’’’’’’’’’’’’’
R
ConstConst*
valueB*  ?*
dtype0*
_output_shapes
:

dense3_sample_weightsPlaceholderWithDefaultConst*
dtype0*#
_output_shapes
:’’’’’’’’’*
shape:’’’’’’’’’
v
total/Initializer/zerosConst*
valueB
 *    *
_class

loc:@total*
dtype0*
_output_shapes
: 
x
totalVarHandleOp*
_class

loc:@total*
dtype0*
_output_shapes
: *
shape: *
shared_nametotal
[
&total/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal*
_output_shapes
: 
g
total/AssignAssignVariableOptotaltotal/Initializer/zeros*
_class

loc:@total*
dtype0
q
total/Read/ReadVariableOpReadVariableOptotal*
_class

loc:@total*
dtype0*
_output_shapes
: 
v
count/Initializer/zerosConst*
valueB
 *    *
_class

loc:@count*
dtype0*
_output_shapes
: 
x
countVarHandleOp*
_class

loc:@count*
dtype0*
_output_shapes
: *
shape: *
shared_namecount
[
&count/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount*
_output_shapes
: 
g
count/AssignAssignVariableOpcountcount/Initializer/zeros*
dtype0*
_class

loc:@count
q
count/Read/ReadVariableOpReadVariableOpcount*
_class

loc:@count*
dtype0*
_output_shapes
: 
\
metrics/accuracy/Cast/xConst*
valueB
 *   ?*
dtype0*
_output_shapes
: 
~
metrics/accuracy/GreaterGreaterdense3/Sigmoidmetrics/accuracy/Cast/x*
T0*'
_output_shapes
:’’’’’’’’’
z
metrics/accuracy/Cast_1Castmetrics/accuracy/Greater*'
_output_shapes
:’’’’’’’’’*

DstT0*

SrcT0


metrics/accuracy/EqualEqualdense3_targetmetrics/accuracy/Cast_1*
T0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’

metrics/accuracy/Cast_2Castmetrics/accuracy/Equal*

SrcT0
*0
_output_shapes
:’’’’’’’’’’’’’’’’’’*

DstT0
r
'metrics/accuracy/Mean/reduction_indicesConst*
valueB :
’’’’’’’’’*
dtype0*
_output_shapes
: 

metrics/accuracy/MeanMeanmetrics/accuracy/Cast_2'metrics/accuracy/Mean/reduction_indices*#
_output_shapes
:’’’’’’’’’*
T0
`
metrics/accuracy/ConstConst*
valueB: *
dtype0*
_output_shapes
:
k
metrics/accuracy/SumSummetrics/accuracy/Meanmetrics/accuracy/Const*
_output_shapes
: *
T0
e
$metrics/accuracy/AssignAddVariableOpAssignAddVariableOptotalmetrics/accuracy/Sum*
dtype0

metrics/accuracy/ReadVariableOpReadVariableOptotal%^metrics/accuracy/AssignAddVariableOp^metrics/accuracy/Sum*
dtype0*
_output_shapes
: 
U
metrics/accuracy/SizeSizemetrics/accuracy/Mean*
T0*
_output_shapes
: 
f
metrics/accuracy/Cast_3Castmetrics/accuracy/Size*

SrcT0*
_output_shapes
: *

DstT0

&metrics/accuracy/AssignAddVariableOp_1AssignAddVariableOpcountmetrics/accuracy/Cast_3%^metrics/accuracy/AssignAddVariableOp*
dtype0
Æ
!metrics/accuracy/ReadVariableOp_1ReadVariableOpcount%^metrics/accuracy/AssignAddVariableOp'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

*metrics/accuracy/div_no_nan/ReadVariableOpReadVariableOptotal'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

,metrics/accuracy/div_no_nan/ReadVariableOp_1ReadVariableOpcount'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
¢
metrics/accuracy/div_no_nanDivNoNan*metrics/accuracy/div_no_nan/ReadVariableOp,metrics/accuracy/div_no_nan/ReadVariableOp_1*
T0*
_output_shapes
: 
c
metrics/accuracy/IdentityIdentitymetrics/accuracy/div_no_nan*
T0*
_output_shapes
: 
^
metrics/accuracy/Cast_4/xConst*
valueB
 *   ?*
dtype0*
_output_shapes
: 

metrics/accuracy/Greater_1Greaterdense3/Sigmoidmetrics/accuracy/Cast_4/x*
T0*'
_output_shapes
:’’’’’’’’’
|
metrics/accuracy/Cast_5Castmetrics/accuracy/Greater_1*

SrcT0
*'
_output_shapes
:’’’’’’’’’*

DstT0

metrics/accuracy/Equal_1Equaldense3_targetmetrics/accuracy/Cast_5*0
_output_shapes
:’’’’’’’’’’’’’’’’’’*
T0

metrics/accuracy/Cast_6Castmetrics/accuracy/Equal_1*

SrcT0
*0
_output_shapes
:’’’’’’’’’’’’’’’’’’*

DstT0
t
)metrics/accuracy/Mean_1/reduction_indicesConst*
dtype0*
_output_shapes
: *
valueB :
’’’’’’’’’

metrics/accuracy/Mean_1Meanmetrics/accuracy/Cast_6)metrics/accuracy/Mean_1/reduction_indices*
T0*#
_output_shapes
:’’’’’’’’’
b
metrics/accuracy/Const_1Const*
valueB: *
dtype0*
_output_shapes
:
s
metrics/accuracy/Mean_2Meanmetrics/accuracy/Mean_1metrics/accuracy/Const_1*
T0*
_output_shapes
: 
¤
5loss/dense3_loss/mean_squared_error/SquaredDifferenceSquaredDifferencedense3/Sigmoiddense3_target*
T0*0
_output_shapes
:’’’’’’’’’’’’’’’’’’

:loss/dense3_loss/mean_squared_error/Mean/reduction_indicesConst*
valueB :
’’’’’’’’’*
dtype0*
_output_shapes
: 
Ń
(loss/dense3_loss/mean_squared_error/MeanMean5loss/dense3_loss/mean_squared_error/SquaredDifference:loss/dense3_loss/mean_squared_error/Mean/reduction_indices*
T0*#
_output_shapes
:’’’’’’’’’
«
floss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapeShapedense3_sample_weights*
T0*
_output_shapes
:
§
eloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rankConst*
dtype0*
_output_shapes
: *
value	B :
½
eloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapeShape(loss/dense3_loss/mean_squared_error/Mean*
T0*
_output_shapes
:
¦
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rankConst*
dtype0*
_output_shapes
: *
value	B :
¦
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar/xConst*
value	B : *
dtype0*
_output_shapes
: 
Ł
bloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarEqualdloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar/xeloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank*
T0*
_output_shapes
: 
ć
nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/SwitchSwitchbloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarbloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar*
T0
*
_output_shapes
: : 

ploss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_tIdentityploss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch:1*
T0
*
_output_shapes
: 

ploss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_fIdentitynloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch*
T0
*
_output_shapes
: 

oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_idIdentitybloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar*
T0
*
_output_shapes
: 
é
ploss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1Switchbloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalaroloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
_output_shapes
: : *
T0
*u
_classk
igloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar
ė
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankEqualloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switchloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1*
T0*
_output_shapes
: 

loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/SwitchSwitchdloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rankoloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
T0*w
_classm
kiloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank*
_output_shapes
: : 

loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1Switcheloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rankoloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
_output_shapes
: : *
T0*x
_classn
ljloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank
Ų
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/SwitchSwitchloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
T0
*
_output_shapes
: : 
Å
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_tIdentityloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:1*
T0
*
_output_shapes
: 
Ć
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_fIdentityloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch*
T0
*
_output_shapes
: 
Č
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_idIdentityloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
T0
*
_output_shapes
: 
ū
”loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dimConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
valueB :
’’’’’’’’’*
dtype0*
_output_shapes
: 
¤
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims
ExpandDimsØloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1”loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim*
T0*
_output_shapes

:
¬
¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/SwitchSwitcheloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapeoloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
T0*x
_classn
ljloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::

¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1Switch¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switchloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0*x
_classn
ljloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::

¢loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/ShapeConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
dtype0*
_output_shapes
:*
valueB"      
ó
¢loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/ConstConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
dtype0*
_output_shapes
: *
value	B :

loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_likeFill¢loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Shape¢loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const*
_output_shapes

:*
T0
ļ
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axisConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
dtype0*
_output_shapes
: *
value	B :
ø
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concatConcatV2loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDimsloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_likeloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis*
T0*
N*
_output_shapes

:
ż
£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dimConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
valueB :
’’’’’’’’’*
dtype0*
_output_shapes
: 
Ŗ
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1
ExpandDimsŖloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim*
T0*
_output_shapes

:
°
¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/SwitchSwitchfloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapeoloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
T0*y
_classo
mkloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::

Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1Switch¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switchloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0*y
_classo
mkloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::
å
«loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperationDenseToDenseSetOperationloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat*
T0*<
_output_shapes*
(:’’’’’’’’’:’’’’’’’’’:*
set_operationa-b
ż
£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dimsSize­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1*
T0*
_output_shapes
: 
å
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/xConst^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
value	B : *
dtype0*
_output_shapes
: 
ś
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dimsEqualloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims*
T0*
_output_shapes
: 
ü
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1Switchloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0
*¤
_class
loc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
_output_shapes
: : 
ß
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/MergeMergeloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims*
T0
*
N*
_output_shapes
: : 
 
mloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/MergeMergeloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Mergerloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1*
N*
_output_shapes
: : *
T0

Ę
^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/ConstConst*
dtype0*
_output_shapes
: *8
value/B- B'weights can not be broadcast to values.
Æ
`loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_1Const*
valueB Bweights.shape=*
dtype0*
_output_shapes
: 
ø
`loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_2Const*(
valueB Bdense3_sample_weights:0*
dtype0*
_output_shapes
: 
®
`loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_3Const*
dtype0*
_output_shapes
: *
valueB Bvalues.shape=
Ė
`loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_4Const*;
value2B0 B*loss/dense3_loss/mean_squared_error/Mean:0*
dtype0*
_output_shapes
: 
«
`loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_5Const*
dtype0*
_output_shapes
: *
valueB B
is_scalar=
ö
kloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/SwitchSwitchmloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Mergemloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge*
_output_shapes
: : *
T0


mloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_tIdentitymloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Switch:1*
T0
*
_output_shapes
: 

mloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_fIdentitykloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Switch*
T0
*
_output_shapes
: 

lloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_idIdentitymloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge*
T0
*
_output_shapes
: 
į
iloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/NoOpNoOpn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t

wloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependencyIdentitymloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_tj^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/NoOp*
T0
*
_classv
trloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t*
_output_shapes
: 
Ź
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*8
value/B- B'weights can not be broadcast to values.*
dtype0*
_output_shapes
: 
±
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
valueB Bweights.shape=*
dtype0*
_output_shapes
: 
ŗ
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*(
valueB Bdense3_sample_weights:0*
dtype0*
_output_shapes
: 
°
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
valueB Bvalues.shape=*
dtype0*
_output_shapes
: 
Ķ
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*;
value2B0 B*loss/dense3_loss/mean_squared_error/Mean:0*
dtype0*
_output_shapes
: 
­
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7Constn^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
valueB B
is_scalar=*
dtype0*
_output_shapes
: 


kloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/AssertAssertrloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switchrloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3*
T
2	

’
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/SwitchSwitchmloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Mergelloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
_output_shapes
: : *
T0
*
_classv
trloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge
ś
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1Switchfloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapelloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0*y
_classo
mkloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::
ų
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2Switcheloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapelloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0*x
_classn
ljloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::
ź
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3Switchbloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarlloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0
*u
_classk
igloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar*
_output_shapes
: : 

yloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1Identitymloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_fl^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert*
_output_shapes
: *
T0
*
_classv
trloc:@loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f

jloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/MergeMergeyloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1wloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency*
T0
*
N*
_output_shapes
: : 

Sloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ShapeShape(loss/dense3_loss/mean_squared_error/Meank^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Merge*
T0*
_output_shapes
:

Sloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ConstConstk^loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Merge*
valueB
 *  ?*
dtype0*
_output_shapes
: 
­
Mloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_likeFillSloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ShapeSloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/Const*
T0*#
_output_shapes
:’’’’’’’’’
Ž
Closs/dense3_loss/mean_squared_error/weighted_loss/broadcast_weightsMuldense3_sample_weightsMloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like*#
_output_shapes
:’’’’’’’’’*
T0
Ł
5loss/dense3_loss/mean_squared_error/weighted_loss/MulMul(loss/dense3_loss/mean_squared_error/MeanCloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights*#
_output_shapes
:’’’’’’’’’*
T0
`
loss/dense3_loss/ConstConst*
valueB: *
dtype0*
_output_shapes
:

loss/dense3_loss/SumSum5loss/dense3_loss/mean_squared_error/weighted_loss/Mulloss/dense3_loss/Const*
_output_shapes
: *
T0
}
loss/dense3_loss/num_elementsSize5loss/dense3_loss/mean_squared_error/weighted_loss/Mul*
T0*
_output_shapes
: 
y
"loss/dense3_loss/num_elements/CastCastloss/dense3_loss/num_elements*
_output_shapes
: *

DstT0*

SrcT0
[
loss/dense3_loss/mul/xConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 
x
loss/dense3_loss/mulMulloss/dense3_loss/mul/x"loss/dense3_loss/num_elements/Cast*
T0*
_output_shapes
: 
[
loss/dense3_loss/Const_1Const*
valueB *
dtype0*
_output_shapes
: 
n
loss/dense3_loss/Sum_1Sumloss/dense3_loss/Sumloss/dense3_loss/Const_1*
T0*
_output_shapes
: 
q
loss/dense3_loss/valueDivNoNanloss/dense3_loss/Sum_1loss/dense3_loss/mul*
T0*
_output_shapes
: 
O

loss/mul/xConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 
T
loss/mulMul
loss/mul/xloss/dense3_loss/value*
_output_shapes
: *
T0
q
iter/Initializer/zerosConst*
dtype0	*
_output_shapes
: *
value	B	 R *
_class
	loc:@iter

iterVarHandleOp"/device:CPU:0*
shared_nameiter*
_class
	loc:@iter*
dtype0	*
_output_shapes
: *
shape: 
h
%iter/IsInitialized/VarIsInitializedOpVarIsInitializedOpiter"/device:CPU:0*
_output_shapes
: 
r
iter/AssignAssignVariableOpiteriter/Initializer/zeros"/device:CPU:0*
dtype0	*
_class
	loc:@iter
}
iter/Read/ReadVariableOpReadVariableOpiter"/device:CPU:0*
_class
	loc:@iter*
dtype0	*
_output_shapes
: 

'learning_rate/Initializer/initial_valueConst*
valueB
 *o:* 
_class
loc:@learning_rate*
dtype0*
_output_shapes
: 

learning_rateVarHandleOp*
shared_namelearning_rate* 
_class
loc:@learning_rate*
dtype0*
_output_shapes
: *
shape: 
k
.learning_rate/IsInitialized/VarIsInitializedOpVarIsInitializedOplearning_rate*
_output_shapes
: 

learning_rate/AssignAssignVariableOplearning_rate'learning_rate/Initializer/initial_value* 
_class
loc:@learning_rate*
dtype0

!learning_rate/Read/ReadVariableOpReadVariableOplearning_rate* 
_class
loc:@learning_rate*
dtype0*
_output_shapes
: 
~
decay/Initializer/initial_valueConst*
dtype0*
_output_shapes
: *
valueB
 *    *
_class

loc:@decay
x
decayVarHandleOp*
shape: *
shared_namedecay*
_class

loc:@decay*
dtype0*
_output_shapes
: 
[
&decay/IsInitialized/VarIsInitializedOpVarIsInitializedOpdecay*
_output_shapes
: 
o
decay/AssignAssignVariableOpdecaydecay/Initializer/initial_value*
_class

loc:@decay*
dtype0
q
decay/Read/ReadVariableOpReadVariableOpdecay*
_class

loc:@decay*
dtype0*
_output_shapes
: 

 beta_1/Initializer/initial_valueConst*
valueB
 *fff?*
_class
loc:@beta_1*
dtype0*
_output_shapes
: 
{
beta_1VarHandleOp*
shape: *
shared_namebeta_1*
_class
loc:@beta_1*
dtype0*
_output_shapes
: 
]
'beta_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpbeta_1*
_output_shapes
: 
s
beta_1/AssignAssignVariableOpbeta_1 beta_1/Initializer/initial_value*
_class
loc:@beta_1*
dtype0
t
beta_1/Read/ReadVariableOpReadVariableOpbeta_1*
_class
loc:@beta_1*
dtype0*
_output_shapes
: 

 beta_2/Initializer/initial_valueConst*
valueB
 *w¾?*
_class
loc:@beta_2*
dtype0*
_output_shapes
: 
{
beta_2VarHandleOp*
shape: *
shared_namebeta_2*
_class
loc:@beta_2*
dtype0*
_output_shapes
: 
]
'beta_2/IsInitialized/VarIsInitializedOpVarIsInitializedOpbeta_2*
_output_shapes
: 
s
beta_2/AssignAssignVariableOpbeta_2 beta_2/Initializer/initial_value*
_class
loc:@beta_2*
dtype0
t
beta_2/Read/ReadVariableOpReadVariableOpbeta_2*
_class
loc:@beta_2*
dtype0*
_output_shapes
: 

!epsilon/Initializer/initial_valueConst*
valueB
 *æÖ3*
_class
loc:@epsilon*
dtype0*
_output_shapes
: 
~
epsilonVarHandleOp*
shared_name	epsilon*
_class
loc:@epsilon*
dtype0*
_output_shapes
: *
shape: 
_
(epsilon/IsInitialized/VarIsInitializedOpVarIsInitializedOpepsilon*
_output_shapes
: 
w
epsilon/AssignAssignVariableOpepsilon!epsilon/Initializer/initial_value*
_class
loc:@epsilon*
dtype0
w
epsilon/Read/ReadVariableOpReadVariableOpepsilon*
_class
loc:@epsilon*
dtype0*
_output_shapes
: 
B
evaluation/group_depsNoOp	^loss/mul^metrics/accuracy/Mean_2
W
Const_1Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_2Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_3Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_4Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_5Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_6Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_7Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
R
VarIsInitializedOpVarIsInitializedOpin1_dense/kernel*
_output_shapes
: 
T
VarIsInitializedOp_1VarIsInitializedOpin2_dense/kernel*
_output_shapes
: 
J
VarIsInitializedOp_2VarIsInitializedOpbeta_2*
_output_shapes
: 
I
VarIsInitializedOp_3VarIsInitializedOptotal*
_output_shapes
: 
R
VarIsInitializedOp_4VarIsInitializedOpin2_dense/bias*
_output_shapes
: 
Q
VarIsInitializedOp_5VarIsInitializedOplearning_rate*
_output_shapes
: 
I
VarIsInitializedOp_6VarIsInitializedOpdecay*
_output_shapes
: 
H
VarIsInitializedOp_7VarIsInitializedOpiter*
_output_shapes
: 
I
VarIsInitializedOp_8VarIsInitializedOpcount*
_output_shapes
: 
K
VarIsInitializedOp_9VarIsInitializedOpepsilon*
_output_shapes
: 
K
VarIsInitializedOp_10VarIsInitializedOpbeta_1*
_output_shapes
: 
R
VarIsInitializedOp_11VarIsInitializedOpdense3/kernel*
_output_shapes
: 
S
VarIsInitializedOp_12VarIsInitializedOpin1_dense/bias*
_output_shapes
: 
P
VarIsInitializedOp_13VarIsInitializedOpdense3/bias*
_output_shapes
: 

	init/NoOpNoOp^beta_1/Assign^beta_2/Assign^count/Assign^decay/Assign^dense3/bias/Assign^dense3/kernel/Assign^epsilon/Assign^in1_dense/bias/Assign^in1_dense/kernel/Assign^in2_dense/bias/Assign^in2_dense/kernel/Assign^learning_rate/Assign^total/Assign
0
init/NoOp_1NoOp^iter/Assign"/device:CPU:0
&
initNoOp
^init/NoOp^init/NoOp_1
W
Const_8Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
\
Const_9Const"/device:CPU:0*
valueB Bmodel*
dtype0*
_output_shapes
: 
X
Const_10Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_11Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_12Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_13Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_14Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_15Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
X
Const_16Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_17Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
¤
RestoreV2/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*K
valueBB@B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
r
RestoreV2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

	RestoreV2	RestoreV2Const_9RestoreV2/tensor_namesRestoreV2/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
B
IdentityIdentity	RestoreV2*
T0*
_output_shapes
:
M
AssignVariableOpAssignVariableOpin1_dense/kernelIdentity*
dtype0
¤
RestoreV2_1/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_1	RestoreV2Const_9RestoreV2_1/tensor_namesRestoreV2_1/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_1IdentityRestoreV2_1*
T0*
_output_shapes
:
O
AssignVariableOp_1AssignVariableOpin1_dense/bias
Identity_1*
dtype0
¦
RestoreV2_2/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_2	RestoreV2Const_9RestoreV2_2/tensor_namesRestoreV2_2/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_2IdentityRestoreV2_2*
T0*
_output_shapes
:
Q
AssignVariableOp_2AssignVariableOpin2_dense/kernel
Identity_2*
dtype0
¤
RestoreV2_3/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_3/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_3	RestoreV2Const_9RestoreV2_3/tensor_namesRestoreV2_3/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_3IdentityRestoreV2_3*
T0*
_output_shapes
:
O
AssignVariableOp_3AssignVariableOpin2_dense/bias
Identity_3*
dtype0
¦
RestoreV2_4/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_4/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_4	RestoreV2Const_9RestoreV2_4/tensor_namesRestoreV2_4/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_4IdentityRestoreV2_4*
T0*
_output_shapes
:
N
AssignVariableOp_4AssignVariableOpdense3/kernel
Identity_4*
dtype0
¤
RestoreV2_5/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_5/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_5	RestoreV2Const_9RestoreV2_5/tensor_namesRestoreV2_5/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_5IdentityRestoreV2_5*
T0*
_output_shapes
:
L
AssignVariableOp_5AssignVariableOpdense3/bias
Identity_5*
dtype0

RestoreV2_6/tensor_namesConst"/device:CPU:0*>
value5B3B)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_6/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_6	RestoreV2Const_9RestoreV2_6/tensor_namesRestoreV2_6/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2	
U

Identity_6IdentityRestoreV2_6"/device:CPU:0*
T0	*
_output_shapes
:
T
AssignVariableOp_6AssignVariableOpiter
Identity_6"/device:CPU:0*
dtype0	
¢
RestoreV2_7/tensor_namesConst"/device:CPU:0*G
value>B<B2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_7/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_7	RestoreV2Const_9RestoreV2_7/tensor_namesRestoreV2_7/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_7IdentityRestoreV2_7*
T0*
_output_shapes
:
N
AssignVariableOp_7AssignVariableOplearning_rate
Identity_7*
dtype0

RestoreV2_8/tensor_namesConst"/device:CPU:0*?
value6B4B*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_8/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_8	RestoreV2Const_9RestoreV2_8/tensor_namesRestoreV2_8/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_8IdentityRestoreV2_8*
_output_shapes
:*
T0
F
AssignVariableOp_8AssignVariableOpdecay
Identity_8*
dtype0

RestoreV2_9/tensor_namesConst"/device:CPU:0*@
value7B5B+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_9/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_9	RestoreV2Const_9RestoreV2_9/tensor_namesRestoreV2_9/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_9IdentityRestoreV2_9*
T0*
_output_shapes
:
G
AssignVariableOp_9AssignVariableOpbeta_1
Identity_9*
dtype0

RestoreV2_10/tensor_namesConst"/device:CPU:0*@
value7B5B+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_10/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_10	RestoreV2Const_9RestoreV2_10/tensor_namesRestoreV2_10/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
H
Identity_10IdentityRestoreV2_10*
T0*
_output_shapes
:
I
AssignVariableOp_10AssignVariableOpbeta_2Identity_10*
dtype0

RestoreV2_11/tensor_namesConst"/device:CPU:0*A
value8B6B,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
u
RestoreV2_11/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_11	RestoreV2Const_9RestoreV2_11/tensor_namesRestoreV2_11/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
H
Identity_11IdentityRestoreV2_11*
T0*
_output_shapes
:
J
AssignVariableOp_11AssignVariableOpepsilonIdentity_11*
dtype0
z
total_1/Initializer/zerosConst*
valueB
 *    *
_class
loc:@total_1*
dtype0*
_output_shapes
: 
~
total_1VarHandleOp*
shared_name	total_1*
_class
loc:@total_1*
dtype0*
_output_shapes
: *
shape: 
_
(total_1/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal_1*
_output_shapes
: 
o
total_1/AssignAssignVariableOptotal_1total_1/Initializer/zeros*
_class
loc:@total_1*
dtype0
w
total_1/Read/ReadVariableOpReadVariableOptotal_1*
dtype0*
_output_shapes
: *
_class
loc:@total_1
z
count_1/Initializer/zerosConst*
valueB
 *    *
_class
loc:@count_1*
dtype0*
_output_shapes
: 
~
count_1VarHandleOp*
_class
loc:@count_1*
dtype0*
_output_shapes
: *
shape: *
shared_name	count_1
_
(count_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount_1*
_output_shapes
: 
o
count_1/AssignAssignVariableOpcount_1count_1/Initializer/zeros*
_class
loc:@count_1*
dtype0
w
count_1/Read/ReadVariableOpReadVariableOpcount_1*
dtype0*
_output_shapes
: *
_class
loc:@count_1
K
Const_18Const*
valueB *
dtype0*
_output_shapes
: 
N
SumSummetrics/accuracy/Mean_2Const_18*
T0*
_output_shapes
: 
E
AssignAddVariableOpAssignAddVariableOptotal_1Sum*
dtype0
j
ReadVariableOpReadVariableOptotal_1^AssignAddVariableOp^Sum*
dtype0*
_output_shapes
: 
F
SizeConst*
value	B :*
dtype0*
_output_shapes
: 
B
CastCastSize*
_output_shapes
: *

DstT0*

SrcT0
^
AssignAddVariableOp_1AssignAddVariableOpcount_1Cast^AssignAddVariableOp*
dtype0
~
ReadVariableOp_1ReadVariableOpcount_1^AssignAddVariableOp^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
q
div_no_nan/ReadVariableOpReadVariableOptotal_1^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
s
div_no_nan/ReadVariableOp_1ReadVariableOpcount_1^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
o

div_no_nanDivNoNandiv_no_nan/ReadVariableOpdiv_no_nan/ReadVariableOp_1*
T0*
_output_shapes
: 
D
Identity_12Identity
div_no_nan*
T0*
_output_shapes
: 
[
div_no_nan_1/ReadVariableOpReadVariableOptotal_1*
dtype0*
_output_shapes
: 
]
div_no_nan_1/ReadVariableOp_1ReadVariableOpcount_1*
dtype0*
_output_shapes
: 
u
div_no_nan_1DivNoNandiv_no_nan_1/ReadVariableOpdiv_no_nan_1/ReadVariableOp_1*
T0*
_output_shapes
: 
F
Identity_13Identitydiv_no_nan_1*
T0*
_output_shapes
: 
l
metric_op_wrapperConst^AssignAddVariableOp_1*
dtype0*
_output_shapes
: *
valueB 
Y
save/filename/inputConst*
valueB Bmodel*
dtype0*
_output_shapes
: 
n
save/filenamePlaceholderWithDefaultsave/filename/input*
shape: *
dtype0*
_output_shapes
: 
e

save/ConstPlaceholderWithDefaultsave/filename*
dtype0*
_output_shapes
: *
shape: 

save/Const_1Const*
dtype0*
_output_shapes
: *Ž
valueŌBŃ BŹ{"class_name": "Model", "config": {"input_layers": [["in1", 0, 0], ["in2", 0, 0]], "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 3], "dtype": "float32", "name": "in1", "sparse": false}, "inbound_nodes": [], "name": "in1"}, {"class_name": "InputLayer", "config": {"batch_input_shape": [null, 2], "dtype": "float32", "name": "in2", "sparse": false}, "inbound_nodes": [], "name": "in2"}, {"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in1_dense", "trainable": true, "units": 5, "use_bias": true}, "inbound_nodes": [["in1", 0, 0, {}]], "name": "in1_dense"}, {"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in2_dense", "trainable": true, "units": 10, "use_bias": true}, "inbound_nodes": [["in2", 0, 0, {}]], "name": "in2_dense"}, {"class_name": "Concatenate", "config": {"axis": -1, "dtype": "float32", "name": "merge", "trainable": true}, "inbound_nodes": [[["in1_dense", 0, 0, {}], ["in2_dense", 0, 0, {}]]], "name": "merge"}, {"class_name": "Dense", "config": {"activation": "sigmoid", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense3", "trainable": true, "units": 1, "use_bias": true}, "inbound_nodes": [["merge", 0, 0, {}]], "name": "dense3"}], "name": "model", "output_layers": ["dense3", 0, 0]}}
Ģ
save/Const_2Const*
valueB B|{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 3], "dtype": "float32", "name": "in1", "sparse": false}}*
dtype0*
_output_shapes
: 
Ģ
save/Const_3Const*
valueB B|{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 2], "dtype": "float32", "name": "in2", "sparse": false}}*
dtype0*
_output_shapes
: 
ŗ
save/Const_4Const*~
valueuBs Bm{"class_name": "Concatenate", "config": {"axis": -1, "dtype": "float32", "name": "merge", "trainable": true}}*
dtype0*
_output_shapes
: 
ų
save/Const_5Const*»
value±B® B§{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in1_dense", "trainable": true, "units": 5, "use_bias": true}}*
dtype0*
_output_shapes
: 
ł
save/Const_6Const*¼
value²BÆ BØ{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in2_dense", "trainable": true, "units": 10, "use_bias": true}}*
dtype0*
_output_shapes
: 
ų
save/Const_7Const*»
value±B® B§{"class_name": "Dense", "config": {"activation": "sigmoid", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense3", "trainable": true, "units": 1, "use_bias": true}}*
dtype0*
_output_shapes
: 
N
save/VarIsInitializedOpVarIsInitializedOptotal_1*
_output_shapes
: 
P
save/VarIsInitializedOp_1VarIsInitializedOpcount_1*
_output_shapes
: 
3
	save/initNoOp^count_1/Assign^total_1/Assign
Ŗ
save/Const_8Const*ķ
valuećBą BŁ{"class_name": "Adam", "config": {"amsgrad": false, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "decay": 0.0, "epsilon": 1.0000000116860974e-07, "learning_rate": 0.0010000000474974513, "name": "Adam"}}*
dtype0*
_output_shapes
: 
”
save/SaveV2/tensor_namesConst*Ō
valueŹBĒB/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-4/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-2/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB(optimizer/.ATTRIBUTES/OBJECT_CONFIG_JSONB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:

save/SaveV2/shape_and_slicesConst*;
value2B0B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:
ņ
save/SaveV2SaveV2
save/Constsave/SaveV2/tensor_namessave/SaveV2/shape_and_slicessave/Const_1save/Const_2save/Const_3save/Const_4save/Const_5"in1_dense/bias/Read/ReadVariableOp$in1_dense/kernel/Read/ReadVariableOpsave/Const_6"in2_dense/bias/Read/ReadVariableOp$in2_dense/kernel/Read/ReadVariableOpsave/Const_7dense3/bias/Read/ReadVariableOp!dense3/kernel/Read/ReadVariableOpsave/Const_8beta_1/Read/ReadVariableOpbeta_2/Read/ReadVariableOpdecay/Read/ReadVariableOpepsilon/Read/ReadVariableOpiter/Read/ReadVariableOp!learning_rate/Read/ReadVariableOp*"
dtypes
2	
}
save/control_dependencyIdentity
save/Const^save/SaveV2*
T0*
_class
loc:@save/Const*
_output_shapes
: 
³
save/RestoreV2/tensor_namesConst"/device:CPU:0*Ō
valueŹBĒB/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-4/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-2/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB(optimizer/.ATTRIBUTES/OBJECT_CONFIG_JSONB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:

save/RestoreV2/shape_and_slicesConst"/device:CPU:0*;
value2B0B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:
ž
save/RestoreV2	RestoreV2
save/Constsave/RestoreV2/tensor_namessave/RestoreV2/shape_and_slices"/device:CPU:0*d
_output_shapesR
P::::::::::::::::::::*"
dtypes
2	

	save/NoOpNoOp

save/NoOp_1NoOp

save/NoOp_2NoOp

save/NoOp_3NoOp

save/NoOp_4NoOp
N
save/IdentityIdentitysave/RestoreV2:5*
T0*
_output_shapes
:
U
save/AssignVariableOpAssignVariableOpin1_dense/biassave/Identity*
dtype0
P
save/Identity_1Identitysave/RestoreV2:6*
T0*
_output_shapes
:
[
save/AssignVariableOp_1AssignVariableOpin1_dense/kernelsave/Identity_1*
dtype0

save/NoOp_5NoOp
P
save/Identity_2Identitysave/RestoreV2:8*
T0*
_output_shapes
:
Y
save/AssignVariableOp_2AssignVariableOpin2_dense/biassave/Identity_2*
dtype0
P
save/Identity_3Identitysave/RestoreV2:9*
T0*
_output_shapes
:
[
save/AssignVariableOp_3AssignVariableOpin2_dense/kernelsave/Identity_3*
dtype0

save/NoOp_6NoOp
Q
save/Identity_4Identitysave/RestoreV2:11*
T0*
_output_shapes
:
V
save/AssignVariableOp_4AssignVariableOpdense3/biassave/Identity_4*
dtype0
Q
save/Identity_5Identitysave/RestoreV2:12*
_output_shapes
:*
T0
X
save/AssignVariableOp_5AssignVariableOpdense3/kernelsave/Identity_5*
dtype0

save/NoOp_7NoOp
Q
save/Identity_6Identitysave/RestoreV2:14*
_output_shapes
:*
T0
Q
save/AssignVariableOp_6AssignVariableOpbeta_1save/Identity_6*
dtype0
Q
save/Identity_7Identitysave/RestoreV2:15*
_output_shapes
:*
T0
Q
save/AssignVariableOp_7AssignVariableOpbeta_2save/Identity_7*
dtype0
Q
save/Identity_8Identitysave/RestoreV2:16*
_output_shapes
:*
T0
P
save/AssignVariableOp_8AssignVariableOpdecaysave/Identity_8*
dtype0
Q
save/Identity_9Identitysave/RestoreV2:17*
T0*
_output_shapes
:
R
save/AssignVariableOp_9AssignVariableOpepsilonsave/Identity_9*
dtype0
a
save/Identity_10Identitysave/RestoreV2:18"/device:CPU:0*
T0	*
_output_shapes
:
`
save/AssignVariableOp_10AssignVariableOpitersave/Identity_10"/device:CPU:0*
dtype0	
R
save/Identity_11Identitysave/RestoreV2:19*
_output_shapes
:*
T0
Z
save/AssignVariableOp_11AssignVariableOplearning_ratesave/Identity_11*
dtype0
Ø
save/restore_all/NoOpNoOp^save/AssignVariableOp^save/AssignVariableOp_1^save/AssignVariableOp_11^save/AssignVariableOp_2^save/AssignVariableOp_3^save/AssignVariableOp_4^save/AssignVariableOp_5^save/AssignVariableOp_6^save/AssignVariableOp_7^save/AssignVariableOp_8^save/AssignVariableOp_9
^save/NoOp^save/NoOp_1^save/NoOp_2^save/NoOp_3^save/NoOp_4^save/NoOp_5^save/NoOp_6^save/NoOp_7
I
save/restore_all/NoOp_1NoOp^save/AssignVariableOp_10"/device:CPU:0
J
save/restore_allNoOp^save/restore_all/NoOp^save/restore_all/NoOp_1
0
init_1NoOp^count_1/Assign^total_1/Assign"D
save/Const:0save/control_dependency:0save/restore_all 5 @F8"
trainable_variablesöó

in1_dense/kernel:0in1_dense/kernel/Assign&in1_dense/kernel/Read/ReadVariableOp:0(2-in1_dense/kernel/Initializer/random_uniform:08
w
in1_dense/bias:0in1_dense/bias/Assign$in1_dense/bias/Read/ReadVariableOp:0(2"in1_dense/bias/Initializer/zeros:08

in2_dense/kernel:0in2_dense/kernel/Assign&in2_dense/kernel/Read/ReadVariableOp:0(2-in2_dense/kernel/Initializer/random_uniform:08
w
in2_dense/bias:0in2_dense/bias/Assign$in2_dense/bias/Read/ReadVariableOp:0(2"in2_dense/bias/Initializer/zeros:08
|
dense3/kernel:0dense3/kernel/Assign#dense3/kernel/Read/ReadVariableOp:0(2*dense3/kernel/Initializer/random_uniform:08
k
dense3/bias:0dense3/bias/Assign!dense3/bias/Read/ReadVariableOp:0(2dense3/bias/Initializer/zeros:08"Ķ
local_variables¹¶
Y
	total_1:0total_1/Assigntotal_1/Read/ReadVariableOp:0(2total_1/Initializer/zeros:0
Y
	count_1:0count_1/Assigncount_1/Read/ReadVariableOp:0(2count_1/Initializer/zeros:0"ö
cond_contextäą
Ś

qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/cond_textqloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_t:0 *ū
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:0
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1
qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_t:0ę
qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0Ś
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1
īq
sloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/cond_text_1qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_f:0*Ė5
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Merge:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Merge:1
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:1
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:1
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:0
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:2
¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1
£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims:0
Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
Ŗloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1
„loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim:0
”loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1:0
 loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat:0
„loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims:0
¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const:0
¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Shape:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0
qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
rloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_f:0
floss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank:0
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1:0
floss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch:0
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0ę
qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0qloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:02ł+
ö+
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/cond_textloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0 *Ę(
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:0
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1
­loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:2
¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1
£loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims:0
Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
Ŗloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1
„loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim:0
”loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1:0
 loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat:0
„loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims:0
¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const:0
¤loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Shape:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0Ŗloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1Ö
Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0Ņ
¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0¦loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0Øloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:12Å
Ā
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/cond_text_1loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0*

loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:1
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0¢
loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0

nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/cond_textnloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t:0 *æ
yloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency:0
nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0
oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t:0ą
nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0
÷
ploss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/cond_text_1nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f:0*”
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch:0
vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1:0
vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2:0
vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5:0
tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7:0
{loss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1:0
nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0
oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f:0
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0
oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge:0
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0į
gloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2:0ē
oloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge:0tloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch:0ą
nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0nloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0ā
hloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1:0Ž
dloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0vloss/dense3_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3:0"`
global_stepQO
M
iter:0iter/Assigniter/Read/ReadVariableOp:0(2iter/Initializer/zeros:0"Ź

	variables¼
¹


in1_dense/kernel:0in1_dense/kernel/Assign&in1_dense/kernel/Read/ReadVariableOp:0(2-in1_dense/kernel/Initializer/random_uniform:08
w
in1_dense/bias:0in1_dense/bias/Assign$in1_dense/bias/Read/ReadVariableOp:0(2"in1_dense/bias/Initializer/zeros:08

in2_dense/kernel:0in2_dense/kernel/Assign&in2_dense/kernel/Read/ReadVariableOp:0(2-in2_dense/kernel/Initializer/random_uniform:08
w
in2_dense/bias:0in2_dense/bias/Assign$in2_dense/bias/Read/ReadVariableOp:0(2"in2_dense/bias/Initializer/zeros:08
|
dense3/kernel:0dense3/kernel/Assign#dense3/kernel/Read/ReadVariableOp:0(2*dense3/kernel/Initializer/random_uniform:08
k
dense3/bias:0dense3/bias/Assign!dense3/bias/Read/ReadVariableOp:0(2dense3/bias/Initializer/zeros:08
M
iter:0iter/Assigniter/Read/ReadVariableOp:0(2iter/Initializer/zeros:0
y
learning_rate:0learning_rate/Assign#learning_rate/Read/ReadVariableOp:0(2)learning_rate/Initializer/initial_value:0
Y
decay:0decay/Assigndecay/Read/ReadVariableOp:0(2!decay/Initializer/initial_value:0
]
beta_1:0beta_1/Assignbeta_1/Read/ReadVariableOp:0(2"beta_1/Initializer/initial_value:0
]
beta_2:0beta_2/Assignbeta_2/Read/ReadVariableOp:0(2"beta_2/Initializer/initial_value:0
a
	epsilon:0epsilon/Assignepsilon/Read/ReadVariableOp:0(2#epsilon/Initializer/initial_value:0*@
__saved_model_init_op'%
__saved_model_init_op
init_1*ō
evalė
@
dense3_target/
dense3_target:0’’’’’’’’’’’’’’’’’’
#
in1
in1:0’’’’’’’’’
#
in2
in2:0’’’’’’’’’
loss

loss/mul:0 9
metrics/accuracy/update_op
metric_op_wrapper:0 -
metrics/accuracy/value
Identity_13:0 =
predictions/dense3'
dense3/Sigmoid:0’’’’’’’’’tensorflow/supervised/evalŁ
å±
:
Add
x"T
y"T
z"T"
Ttype:
2	
B
AssignVariableOp
resource
value"dtype"
dtypetype
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
h
ConcatV2
values"T*N
axis"Tidx
output"T"
Nint(0"	
Ttype"
Tidxtype0:
2	
8
Const
output"dtype"
valuetensor"
dtypetype
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
=
Mul
x"T
y"T
z"T"
Ttype:
2	

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
X
PlaceholderWithDefault
input"dtype
output"dtype"
dtypetype"
shapeshape
~
RandomUniform

shape"T
output"dtype"
seedint "
seed2int "
dtypetype:
2"
Ttype:
2	
@
ReadVariableOp
resource
value"dtype"
dtypetype
E
Relu
features"T
activations"T"
Ttype:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
0
Sigmoid
x"T
y"T"
Ttype:

2
:
Sub
x"T
y"T
z"T"
Ttype:
2	
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape
9
VarIsInitializedOp
resource
is_initialized
"serve*2.0.0-alpha02v1.12.0-9492-g2c319fb4158Ż½
f
in1Placeholder*
dtype0*'
_output_shapes
:’’’’’’’’’*
shape:’’’’’’’’’
f
in2Placeholder*
dtype0*'
_output_shapes
:’’’’’’’’’*
shape:’’’’’’’’’
§
1in1_dense/kernel/Initializer/random_uniform/shapeConst*
valueB"      *#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes
:

/in1_dense/kernel/Initializer/random_uniform/minConst*
valueB
 *×³]æ*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes
: 

/in1_dense/kernel/Initializer/random_uniform/maxConst*
valueB
 *×³]?*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes
: 
Ų
9in1_dense/kernel/Initializer/random_uniform/RandomUniformRandomUniform1in1_dense/kernel/Initializer/random_uniform/shape*
T0*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes

:
Ž
/in1_dense/kernel/Initializer/random_uniform/subSub/in1_dense/kernel/Initializer/random_uniform/max/in1_dense/kernel/Initializer/random_uniform/min*
T0*#
_class
loc:@in1_dense/kernel*
_output_shapes
: 
š
/in1_dense/kernel/Initializer/random_uniform/mulMul9in1_dense/kernel/Initializer/random_uniform/RandomUniform/in1_dense/kernel/Initializer/random_uniform/sub*
_output_shapes

:*
T0*#
_class
loc:@in1_dense/kernel
ā
+in1_dense/kernel/Initializer/random_uniformAdd/in1_dense/kernel/Initializer/random_uniform/mul/in1_dense/kernel/Initializer/random_uniform/min*
T0*#
_class
loc:@in1_dense/kernel*
_output_shapes

:
”
in1_dense/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape
:*!
shared_namein1_dense/kernel*#
_class
loc:@in1_dense/kernel
q
1in1_dense/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpin1_dense/kernel*
_output_shapes
: 

in1_dense/kernel/AssignAssignVariableOpin1_dense/kernel+in1_dense/kernel/Initializer/random_uniform*#
_class
loc:@in1_dense/kernel*
dtype0

$in1_dense/kernel/Read/ReadVariableOpReadVariableOpin1_dense/kernel*#
_class
loc:@in1_dense/kernel*
dtype0*
_output_shapes

:

 in1_dense/bias/Initializer/zerosConst*
valueB*    *!
_class
loc:@in1_dense/bias*
dtype0*
_output_shapes
:

in1_dense/biasVarHandleOp*
shape:*
shared_namein1_dense/bias*!
_class
loc:@in1_dense/bias*
dtype0*
_output_shapes
: 
m
/in1_dense/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpin1_dense/bias*
_output_shapes
: 

in1_dense/bias/AssignAssignVariableOpin1_dense/bias in1_dense/bias/Initializer/zeros*!
_class
loc:@in1_dense/bias*
dtype0

"in1_dense/bias/Read/ReadVariableOpReadVariableOpin1_dense/bias*
_output_shapes
:*!
_class
loc:@in1_dense/bias*
dtype0
p
in1_dense/MatMul/ReadVariableOpReadVariableOpin1_dense/kernel*
dtype0*
_output_shapes

:
r
in1_dense/MatMulMatMulin1in1_dense/MatMul/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’
k
 in1_dense/BiasAdd/ReadVariableOpReadVariableOpin1_dense/bias*
dtype0*
_output_shapes
:

in1_dense/BiasAddBiasAddin1_dense/MatMul in1_dense/BiasAdd/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’
[
in1_dense/ReluReluin1_dense/BiasAdd*'
_output_shapes
:’’’’’’’’’*
T0
§
1in2_dense/kernel/Initializer/random_uniform/shapeConst*
valueB"   
   *#
_class
loc:@in2_dense/kernel*
dtype0*
_output_shapes
:

/in2_dense/kernel/Initializer/random_uniform/minConst*
valueB
 *ó5æ*#
_class
loc:@in2_dense/kernel*
dtype0*
_output_shapes
: 

/in2_dense/kernel/Initializer/random_uniform/maxConst*
valueB
 *ó5?*#
_class
loc:@in2_dense/kernel*
dtype0*
_output_shapes
: 
Ų
9in2_dense/kernel/Initializer/random_uniform/RandomUniformRandomUniform1in2_dense/kernel/Initializer/random_uniform/shape*
T0*#
_class
loc:@in2_dense/kernel*
dtype0*
_output_shapes

:

Ž
/in2_dense/kernel/Initializer/random_uniform/subSub/in2_dense/kernel/Initializer/random_uniform/max/in2_dense/kernel/Initializer/random_uniform/min*
_output_shapes
: *
T0*#
_class
loc:@in2_dense/kernel
š
/in2_dense/kernel/Initializer/random_uniform/mulMul9in2_dense/kernel/Initializer/random_uniform/RandomUniform/in2_dense/kernel/Initializer/random_uniform/sub*
T0*#
_class
loc:@in2_dense/kernel*
_output_shapes

:

ā
+in2_dense/kernel/Initializer/random_uniformAdd/in2_dense/kernel/Initializer/random_uniform/mul/in2_dense/kernel/Initializer/random_uniform/min*
_output_shapes

:
*
T0*#
_class
loc:@in2_dense/kernel
”
in2_dense/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape
:
*!
shared_namein2_dense/kernel*#
_class
loc:@in2_dense/kernel
q
1in2_dense/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpin2_dense/kernel*
_output_shapes
: 

in2_dense/kernel/AssignAssignVariableOpin2_dense/kernel+in2_dense/kernel/Initializer/random_uniform*#
_class
loc:@in2_dense/kernel*
dtype0

$in2_dense/kernel/Read/ReadVariableOpReadVariableOpin2_dense/kernel*
dtype0*
_output_shapes

:
*#
_class
loc:@in2_dense/kernel

 in2_dense/bias/Initializer/zerosConst*
valueB
*    *!
_class
loc:@in2_dense/bias*
dtype0*
_output_shapes
:


in2_dense/biasVarHandleOp*
dtype0*
_output_shapes
: *
shape:
*
shared_namein2_dense/bias*!
_class
loc:@in2_dense/bias
m
/in2_dense/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpin2_dense/bias*
_output_shapes
: 

in2_dense/bias/AssignAssignVariableOpin2_dense/bias in2_dense/bias/Initializer/zeros*!
_class
loc:@in2_dense/bias*
dtype0

"in2_dense/bias/Read/ReadVariableOpReadVariableOpin2_dense/bias*!
_class
loc:@in2_dense/bias*
dtype0*
_output_shapes
:

p
in2_dense/MatMul/ReadVariableOpReadVariableOpin2_dense/kernel*
dtype0*
_output_shapes

:

r
in2_dense/MatMulMatMulin2in2_dense/MatMul/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’

k
 in2_dense/BiasAdd/ReadVariableOpReadVariableOpin2_dense/bias*
dtype0*
_output_shapes
:


in2_dense/BiasAddBiasAddin2_dense/MatMul in2_dense/BiasAdd/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’

[
in2_dense/ReluReluin2_dense/BiasAdd*
T0*'
_output_shapes
:’’’’’’’’’

S
merge/concat/axisConst*
value	B :*
dtype0*
_output_shapes
: 

merge/concatConcatV2in1_dense/Reluin2_dense/Relumerge/concat/axis*
T0*
N*'
_output_shapes
:’’’’’’’’’
”
.dense3/kernel/Initializer/random_uniform/shapeConst*
valueB"      * 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes
:

,dense3/kernel/Initializer/random_uniform/minConst*
valueB
 *qÄæ* 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes
: 

,dense3/kernel/Initializer/random_uniform/maxConst*
valueB
 *qÄ?* 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes
: 
Ļ
6dense3/kernel/Initializer/random_uniform/RandomUniformRandomUniform.dense3/kernel/Initializer/random_uniform/shape*
dtype0*
_output_shapes

:*
T0* 
_class
loc:@dense3/kernel
Ņ
,dense3/kernel/Initializer/random_uniform/subSub,dense3/kernel/Initializer/random_uniform/max,dense3/kernel/Initializer/random_uniform/min*
T0* 
_class
loc:@dense3/kernel*
_output_shapes
: 
ä
,dense3/kernel/Initializer/random_uniform/mulMul6dense3/kernel/Initializer/random_uniform/RandomUniform,dense3/kernel/Initializer/random_uniform/sub*
T0* 
_class
loc:@dense3/kernel*
_output_shapes

:
Ö
(dense3/kernel/Initializer/random_uniformAdd,dense3/kernel/Initializer/random_uniform/mul,dense3/kernel/Initializer/random_uniform/min*
T0* 
_class
loc:@dense3/kernel*
_output_shapes

:

dense3/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape
:*
shared_namedense3/kernel* 
_class
loc:@dense3/kernel
k
.dense3/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense3/kernel*
_output_shapes
: 

dense3/kernel/AssignAssignVariableOpdense3/kernel(dense3/kernel/Initializer/random_uniform* 
_class
loc:@dense3/kernel*
dtype0

!dense3/kernel/Read/ReadVariableOpReadVariableOpdense3/kernel* 
_class
loc:@dense3/kernel*
dtype0*
_output_shapes

:

dense3/bias/Initializer/zerosConst*
dtype0*
_output_shapes
:*
valueB*    *
_class
loc:@dense3/bias

dense3/biasVarHandleOp*
dtype0*
_output_shapes
: *
shape:*
shared_namedense3/bias*
_class
loc:@dense3/bias
g
,dense3/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense3/bias*
_output_shapes
: 

dense3/bias/AssignAssignVariableOpdense3/biasdense3/bias/Initializer/zeros*
dtype0*
_class
loc:@dense3/bias

dense3/bias/Read/ReadVariableOpReadVariableOpdense3/bias*
dtype0*
_output_shapes
:*
_class
loc:@dense3/bias
j
dense3/MatMul/ReadVariableOpReadVariableOpdense3/kernel*
dtype0*
_output_shapes

:
u
dense3/MatMulMatMulmerge/concatdense3/MatMul/ReadVariableOp*
T0*'
_output_shapes
:’’’’’’’’’
e
dense3/BiasAdd/ReadVariableOpReadVariableOpdense3/bias*
dtype0*
_output_shapes
:
y
dense3/BiasAddBiasAdddense3/MatMuldense3/BiasAdd/ReadVariableOp*'
_output_shapes
:’’’’’’’’’*
T0
[
dense3/SigmoidSigmoiddense3/BiasAdd*
T0*'
_output_shapes
:’’’’’’’’’
+
predict/group_depsNoOp^dense3/Sigmoid
U
ConstConst"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_1Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_2Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_3Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_4Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_5Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_6Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
\
Const_7Const"/device:CPU:0*
valueB Bmodel*
dtype0*
_output_shapes
: 
W
Const_8Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_9Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_10Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_11Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_12Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_13Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_14Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
¤
RestoreV2/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*K
valueBB@B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
r
RestoreV2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

	RestoreV2	RestoreV2Const_7RestoreV2/tensor_namesRestoreV2/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
B
IdentityIdentity	RestoreV2*
T0*
_output_shapes
:
M
AssignVariableOpAssignVariableOpin1_dense/kernelIdentity*
dtype0
¤
RestoreV2_1/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_1	RestoreV2Const_7RestoreV2_1/tensor_namesRestoreV2_1/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_1IdentityRestoreV2_1*
T0*
_output_shapes
:
O
AssignVariableOp_1AssignVariableOpin1_dense/bias
Identity_1*
dtype0
¦
RestoreV2_2/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_2/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_2	RestoreV2Const_7RestoreV2_2/tensor_namesRestoreV2_2/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_2IdentityRestoreV2_2*
T0*
_output_shapes
:
Q
AssignVariableOp_2AssignVariableOpin2_dense/kernel
Identity_2*
dtype0
¤
RestoreV2_3/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_3/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_3	RestoreV2Const_7RestoreV2_3/tensor_namesRestoreV2_3/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_3IdentityRestoreV2_3*
T0*
_output_shapes
:
O
AssignVariableOp_3AssignVariableOpin2_dense/bias
Identity_3*
dtype0
¦
RestoreV2_4/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_4/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_4	RestoreV2Const_7RestoreV2_4/tensor_namesRestoreV2_4/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_4IdentityRestoreV2_4*
T0*
_output_shapes
:
N
AssignVariableOp_4AssignVariableOpdense3/kernel
Identity_4*
dtype0
¤
RestoreV2_5/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_5/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_5	RestoreV2Const_7RestoreV2_5/tensor_namesRestoreV2_5/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_5IdentityRestoreV2_5*
T0*
_output_shapes
:
L
AssignVariableOp_5AssignVariableOpdense3/bias
Identity_5*
dtype0
M
VarIsInitializedOpVarIsInitializedOpdense3/bias*
_output_shapes
: 
T
VarIsInitializedOp_1VarIsInitializedOpin2_dense/kernel*
_output_shapes
: 
R
VarIsInitializedOp_2VarIsInitializedOpin1_dense/bias*
_output_shapes
: 
R
VarIsInitializedOp_3VarIsInitializedOpin2_dense/bias*
_output_shapes
: 
T
VarIsInitializedOp_4VarIsInitializedOpin1_dense/kernel*
_output_shapes
: 
Q
VarIsInitializedOp_5VarIsInitializedOpdense3/kernel*
_output_shapes
: 

initNoOp^dense3/bias/Assign^dense3/kernel/Assign^in1_dense/bias/Assign^in1_dense/kernel/Assign^in2_dense/bias/Assign^in2_dense/kernel/Assign
Y
save/filename/inputConst*
valueB Bmodel*
dtype0*
_output_shapes
: 
n
save/filenamePlaceholderWithDefaultsave/filename/input*
dtype0*
_output_shapes
: *
shape: 
e

save/ConstPlaceholderWithDefaultsave/filename*
shape: *
dtype0*
_output_shapes
: 

save/Const_1Const*Ž
valueŌBŃ BŹ{"class_name": "Model", "config": {"input_layers": [["in1", 0, 0], ["in2", 0, 0]], "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 3], "dtype": "float32", "name": "in1", "sparse": false}, "inbound_nodes": [], "name": "in1"}, {"class_name": "InputLayer", "config": {"batch_input_shape": [null, 2], "dtype": "float32", "name": "in2", "sparse": false}, "inbound_nodes": [], "name": "in2"}, {"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in1_dense", "trainable": true, "units": 5, "use_bias": true}, "inbound_nodes": [["in1", 0, 0, {}]], "name": "in1_dense"}, {"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in2_dense", "trainable": true, "units": 10, "use_bias": true}, "inbound_nodes": [["in2", 0, 0, {}]], "name": "in2_dense"}, {"class_name": "Concatenate", "config": {"axis": -1, "dtype": "float32", "name": "merge", "trainable": true}, "inbound_nodes": [[["in1_dense", 0, 0, {}], ["in2_dense", 0, 0, {}]]], "name": "merge"}, {"class_name": "Dense", "config": {"activation": "sigmoid", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense3", "trainable": true, "units": 1, "use_bias": true}, "inbound_nodes": [["merge", 0, 0, {}]], "name": "dense3"}], "name": "model", "output_layers": ["dense3", 0, 0]}}*
dtype0*
_output_shapes
: 
Ģ
save/Const_2Const*
valueB B|{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 3], "dtype": "float32", "name": "in1", "sparse": false}}*
dtype0*
_output_shapes
: 
Ģ
save/Const_3Const*
valueB B|{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 2], "dtype": "float32", "name": "in2", "sparse": false}}*
dtype0*
_output_shapes
: 
ŗ
save/Const_4Const*
dtype0*
_output_shapes
: *~
valueuBs Bm{"class_name": "Concatenate", "config": {"axis": -1, "dtype": "float32", "name": "merge", "trainable": true}}
ų
save/Const_5Const*»
value±B® B§{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in1_dense", "trainable": true, "units": 5, "use_bias": true}}*
dtype0*
_output_shapes
: 
ł
save/Const_6Const*
dtype0*
_output_shapes
: *¼
value²BÆ BØ{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "in2_dense", "trainable": true, "units": 10, "use_bias": true}}
ų
save/Const_7Const*»
value±B® B§{"class_name": "Dense", "config": {"activation": "sigmoid", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense3", "trainable": true, "units": 1, "use_bias": true}}*
dtype0*
_output_shapes
: 
ä
save/SaveV2/tensor_namesConst*
valueBB/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-4/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-2/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
}
save/SaveV2/shape_and_slicesConst*-
value$B"B B B B B B B B B B B B B *
dtype0*
_output_shapes
:
°
save/SaveV2SaveV2
save/Constsave/SaveV2/tensor_namessave/SaveV2/shape_and_slicessave/Const_1save/Const_2save/Const_3save/Const_4save/Const_5"in1_dense/bias/Read/ReadVariableOp$in1_dense/kernel/Read/ReadVariableOpsave/Const_6"in2_dense/bias/Read/ReadVariableOp$in2_dense/kernel/Read/ReadVariableOpsave/Const_7dense3/bias/Read/ReadVariableOp!dense3/kernel/Read/ReadVariableOp*
dtypes
2
}
save/control_dependencyIdentity
save/Const^save/SaveV2*
T0*
_class
loc:@save/Const*
_output_shapes
: 
ö
save/RestoreV2/tensor_namesConst"/device:CPU:0*
valueBB/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-4/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-2/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:

save/RestoreV2/shape_and_slicesConst"/device:CPU:0*-
value$B"B B B B B B B B B B B B B *
dtype0*
_output_shapes
:
Ū
save/RestoreV2	RestoreV2
save/Constsave/RestoreV2/tensor_namessave/RestoreV2/shape_and_slices"/device:CPU:0*
dtypes
2*H
_output_shapes6
4:::::::::::::

	save/NoOpNoOp

save/NoOp_1NoOp

save/NoOp_2NoOp

save/NoOp_3NoOp

save/NoOp_4NoOp
N
save/IdentityIdentitysave/RestoreV2:5*
_output_shapes
:*
T0
U
save/AssignVariableOpAssignVariableOpin1_dense/biassave/Identity*
dtype0
P
save/Identity_1Identitysave/RestoreV2:6*
T0*
_output_shapes
:
[
save/AssignVariableOp_1AssignVariableOpin1_dense/kernelsave/Identity_1*
dtype0

save/NoOp_5NoOp
P
save/Identity_2Identitysave/RestoreV2:8*
T0*
_output_shapes
:
Y
save/AssignVariableOp_2AssignVariableOpin2_dense/biassave/Identity_2*
dtype0
P
save/Identity_3Identitysave/RestoreV2:9*
T0*
_output_shapes
:
[
save/AssignVariableOp_3AssignVariableOpin2_dense/kernelsave/Identity_3*
dtype0

save/NoOp_6NoOp
Q
save/Identity_4Identitysave/RestoreV2:11*
T0*
_output_shapes
:
V
save/AssignVariableOp_4AssignVariableOpdense3/biassave/Identity_4*
dtype0
Q
save/Identity_5Identitysave/RestoreV2:12*
T0*
_output_shapes
:
X
save/AssignVariableOp_5AssignVariableOpdense3/kernelsave/Identity_5*
dtype0

save/restore_allNoOp^save/AssignVariableOp^save/AssignVariableOp_1^save/AssignVariableOp_2^save/AssignVariableOp_3^save/AssignVariableOp_4^save/AssignVariableOp_5
^save/NoOp^save/NoOp_1^save/NoOp_2^save/NoOp_3^save/NoOp_4^save/NoOp_5^save/NoOp_6

init_1NoOp"D
save/Const:0save/control_dependency:0save/restore_all 5 @F8"
	variablesöó

in1_dense/kernel:0in1_dense/kernel/Assign&in1_dense/kernel/Read/ReadVariableOp:0(2-in1_dense/kernel/Initializer/random_uniform:08
w
in1_dense/bias:0in1_dense/bias/Assign$in1_dense/bias/Read/ReadVariableOp:0(2"in1_dense/bias/Initializer/zeros:08

in2_dense/kernel:0in2_dense/kernel/Assign&in2_dense/kernel/Read/ReadVariableOp:0(2-in2_dense/kernel/Initializer/random_uniform:08
w
in2_dense/bias:0in2_dense/bias/Assign$in2_dense/bias/Read/ReadVariableOp:0(2"in2_dense/bias/Initializer/zeros:08
|
dense3/kernel:0dense3/kernel/Assign#dense3/kernel/Read/ReadVariableOp:0(2*dense3/kernel/Initializer/random_uniform:08
k
dense3/bias:0dense3/bias/Assign!dense3/bias/Read/ReadVariableOp:0(2dense3/bias/Initializer/zeros:08"
trainable_variablesöó

in1_dense/kernel:0in1_dense/kernel/Assign&in1_dense/kernel/Read/ReadVariableOp:0(2-in1_dense/kernel/Initializer/random_uniform:08
w
in1_dense/bias:0in1_dense/bias/Assign$in1_dense/bias/Read/ReadVariableOp:0(2"in1_dense/bias/Initializer/zeros:08

in2_dense/kernel:0in2_dense/kernel/Assign&in2_dense/kernel/Read/ReadVariableOp:0(2-in2_dense/kernel/Initializer/random_uniform:08
w
in2_dense/bias:0in2_dense/bias/Assign$in2_dense/bias/Read/ReadVariableOp:0(2"in2_dense/bias/Initializer/zeros:08
|
dense3/kernel:0dense3/kernel/Assign#dense3/kernel/Read/ReadVariableOp:0(2*dense3/kernel/Initializer/random_uniform:08
k
dense3/bias:0dense3/bias/Assign!dense3/bias/Read/ReadVariableOp:0(2dense3/bias/Initializer/zeros:08*­
serving_default
#
in1
in1:0’’’’’’’’’
#
in2
in2:0’’’’’’’’’1
dense3'
dense3/Sigmoid:0’’’’’’’’’tensorflow/serving/predict*@
__saved_model_init_op'%
__saved_model_init_op
init_1