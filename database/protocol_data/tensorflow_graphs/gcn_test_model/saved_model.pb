╗Ѕ
ф§
8
Const
output"dtype"
valuetensor"
dtypetype

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetypeѕ
Й
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ѕ
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshapeѕ"serve*2.2.02v2.2.0-rc4-8-g2b96f3662b8ЉН	
ѕ
Conv1Layer/root_kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*'
shared_nameConv1Layer/root_kernel
Ђ
*Conv1Layer/root_kernel/Read/ReadVariableOpReadVariableOpConv1Layer/root_kernel*
_output_shapes

:*
dtype0
v
Conv1Layer/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameConv1Layer/bias
o
#Conv1Layer/bias/Read/ReadVariableOpReadVariableOpConv1Layer/bias*
_output_shapes
:*
dtype0
ђ
OutputLayer/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*#
shared_nameOutputLayer/kernel
y
&OutputLayer/kernel/Read/ReadVariableOpReadVariableOpOutputLayer/kernel*
_output_shapes

:*
dtype0
x
OutputLayer/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*!
shared_nameOutputLayer/bias
q
$OutputLayer/bias/Read/ReadVariableOpReadVariableOpOutputLayer/bias*
_output_shapes
:*
dtype0
і
Conv1Layer/FGN_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*(
shared_nameConv1Layer/FGN_0/kernel
Ѓ
+Conv1Layer/FGN_0/kernel/Read/ReadVariableOpReadVariableOpConv1Layer/FGN_0/kernel*
_output_shapes

:*
dtype0
ѓ
Conv1Layer/FGN_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*&
shared_nameConv1Layer/FGN_0/bias
{
)Conv1Layer/FGN_0/bias/Read/ReadVariableOpReadVariableOpConv1Layer/FGN_0/bias*
_output_shapes
:*
dtype0
і
Conv1Layer/FGN_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*(
shared_nameConv1Layer/FGN_1/kernel
Ѓ
+Conv1Layer/FGN_1/kernel/Read/ReadVariableOpReadVariableOpConv1Layer/FGN_1/kernel*
_output_shapes

:*
dtype0
ѓ
Conv1Layer/FGN_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*&
shared_nameConv1Layer/FGN_1/bias
{
)Conv1Layer/FGN_1/bias/Read/ReadVariableOpReadVariableOpConv1Layer/FGN_1/bias*
_output_shapes
:*
dtype0
ј
Conv1Layer/FGN_out/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:**
shared_nameConv1Layer/FGN_out/kernel
Є
-Conv1Layer/FGN_out/kernel/Read/ReadVariableOpReadVariableOpConv1Layer/FGN_out/kernel*
_output_shapes

:*
dtype0
є
Conv1Layer/FGN_out/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*(
shared_nameConv1Layer/FGN_out/bias

+Conv1Layer/FGN_out/bias/Read/ReadVariableOpReadVariableOpConv1Layer/FGN_out/bias*
_output_shapes
:*
dtype0

NoOpNoOp
░
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*в
valueрBя BО
§
layer-0
layer-1
layer-2
layer_with_weights-0
layer-3
layer_with_weights-1
layer-4
	optimizer
loss
_layers
	trainable_variables

	variables
regularization_losses
	keras_api

signatures
 
 
 
ѕ
kernel_network_layers
root_kernel
bias
trainable_variables
	variables
regularization_losses
	keras_api
h

kernel
bias
trainable_variables
	variables
regularization_losses
	keras_api
 
 
#
0
1
2
3
4
F
0
1
2
3
4
5
6
 7
8
9
F
0
1
2
3
4
5
6
 7
8
9
 
Г
!layer_metrics
	trainable_variables
"layer_regularization_losses

	variables

#layers
$non_trainable_variables
%metrics
regularization_losses
 

&0
'1
(2
ge
VARIABLE_VALUEConv1Layer/root_kernel;layer_with_weights-0/root_kernel/.ATTRIBUTES/VARIABLE_VALUE
YW
VARIABLE_VALUEConv1Layer/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE
8
0
1
2
3
4
5
6
 7
8
0
1
2
3
4
5
6
 7
 
Г
)layer_metrics
trainable_variables
*layer_regularization_losses
	variables

+layers
,non_trainable_variables
-metrics
regularization_losses
^\
VARIABLE_VALUEOutputLayer/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE
ZX
VARIABLE_VALUEOutputLayer/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE

0
1

0
1
 
Г
.layer_metrics
trainable_variables
/layer_regularization_losses
	variables

0layers
1non_trainable_variables
2metrics
regularization_losses
][
VARIABLE_VALUEConv1Layer/FGN_0/kernel0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUE
[Y
VARIABLE_VALUEConv1Layer/FGN_0/bias0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUE
][
VARIABLE_VALUEConv1Layer/FGN_1/kernel0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUE
[Y
VARIABLE_VALUEConv1Layer/FGN_1/bias0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUE
_]
VARIABLE_VALUEConv1Layer/FGN_out/kernel0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUE
][
VARIABLE_VALUEConv1Layer/FGN_out/bias0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUE
 
 
#
0
1
2
3
4
 
 
h

kernel
bias
3trainable_variables
4	variables
5regularization_losses
6	keras_api
h

kernel
bias
7trainable_variables
8	variables
9regularization_losses
:	keras_api
h

kernel
 bias
;trainable_variables
<	variables
=regularization_losses
>	keras_api
 
 

&0
'1
(2
 
 
 
 
 
 
 

0
1

0
1
 
Г
?layer_metrics
3trainable_variables
@layer_regularization_losses
4	variables

Alayers
Bnon_trainable_variables
Cmetrics
5regularization_losses

0
1

0
1
 
Г
Dlayer_metrics
7trainable_variables
Elayer_regularization_losses
8	variables

Flayers
Gnon_trainable_variables
Hmetrics
9regularization_losses

0
 1

0
 1
 
Г
Ilayer_metrics
;trainable_variables
Jlayer_regularization_losses
<	variables

Klayers
Lnon_trainable_variables
Mmetrics
=regularization_losses
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

serving_default_A_inPlaceholder*+
_output_shapes
:         *
dtype0* 
shape:         
Є
serving_default_E_inPlaceholder*/
_output_shapes
:         *
dtype0*$
shape:         

serving_default_X_inPlaceholder*+
_output_shapes
:         *
dtype0* 
shape:         
│
StatefulPartitionedCallStatefulPartitionedCallserving_default_A_inserving_default_E_inserving_default_X_inConv1Layer/FGN_0/kernelConv1Layer/FGN_0/biasConv1Layer/FGN_1/kernelConv1Layer/FGN_1/biasConv1Layer/FGN_out/kernelConv1Layer/FGN_out/biasConv1Layer/root_kernelConv1Layer/biasOutputLayer/kernelOutputLayer/bias*
Tin
2*
Tout
2*+
_output_shapes
:         *,
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8**
f%R#
!__inference_signature_wrapper_891
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
г
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename*Conv1Layer/root_kernel/Read/ReadVariableOp#Conv1Layer/bias/Read/ReadVariableOp&OutputLayer/kernel/Read/ReadVariableOp$OutputLayer/bias/Read/ReadVariableOp+Conv1Layer/FGN_0/kernel/Read/ReadVariableOp)Conv1Layer/FGN_0/bias/Read/ReadVariableOp+Conv1Layer/FGN_1/kernel/Read/ReadVariableOp)Conv1Layer/FGN_1/bias/Read/ReadVariableOp-Conv1Layer/FGN_out/kernel/Read/ReadVariableOp+Conv1Layer/FGN_out/bias/Read/ReadVariableOpConst*
Tin
2*
Tout
2*
_output_shapes
: * 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*&
f!R
__inference__traced_save_1451
▀
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameConv1Layer/root_kernelConv1Layer/biasOutputLayer/kernelOutputLayer/biasConv1Layer/FGN_0/kernelConv1Layer/FGN_0/biasConv1Layer/FGN_1/kernelConv1Layer/FGN_1/biasConv1Layer/FGN_out/kernelConv1Layer/FGN_out/bias*
Tin
2*
Tout
2*
_output_shapes
: * 
_read_only_resource_inputs
 **
config_proto

CPU

GPU 2J 8*)
f$R"
 __inference__traced_restore_1493№Њ	
Ў└
Ј
A__inference_MyModel_layer_call_and_return_conditional_losses_1165
inputs_0
inputs_1
inputs_26
2conv1layer_fgn_0_tensordot_readvariableop_resource4
0conv1layer_fgn_0_biasadd_readvariableop_resource6
2conv1layer_fgn_1_tensordot_readvariableop_resource4
0conv1layer_fgn_1_biasadd_readvariableop_resource8
4conv1layer_fgn_out_tensordot_readvariableop_resource6
2conv1layer_fgn_out_biasadd_readvariableop_resource-
)conv1layer_matmul_readvariableop_resource.
*conv1layer_biasadd_readvariableop_resource1
-outputlayer_tensordot_readvariableop_resource/
+outputlayer_biasadd_readvariableop_resource
identityѕ\
Conv1Layer/ShapeShapeinputs_0*
T0*
_output_shapes
:2
Conv1Layer/ShapeЊ
Conv1Layer/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2 
Conv1Layer/strided_slice/stackЌ
 Conv1Layer/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2"
 Conv1Layer/strided_slice/stack_1ј
 Conv1Layer/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2"
 Conv1Layer/strided_slice/stack_2ц
Conv1Layer/strided_sliceStridedSliceConv1Layer/Shape:output:0'Conv1Layer/strided_slice/stack:output:0)Conv1Layer/strided_slice/stack_1:output:0)Conv1Layer/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
Conv1Layer/strided_slice╔
)Conv1Layer/FGN_0/Tensordot/ReadVariableOpReadVariableOp2conv1layer_fgn_0_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02+
)Conv1Layer/FGN_0/Tensordot/ReadVariableOpї
Conv1Layer/FGN_0/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2!
Conv1Layer/FGN_0/Tensordot/axesЌ
Conv1Layer/FGN_0/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv1Layer/FGN_0/Tensordot/free|
 Conv1Layer/FGN_0/Tensordot/ShapeShapeinputs_2*
T0*
_output_shapes
:2"
 Conv1Layer/FGN_0/Tensordot/Shapeќ
(Conv1Layer/FGN_0/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2*
(Conv1Layer/FGN_0/Tensordot/GatherV2/axisд
#Conv1Layer/FGN_0/Tensordot/GatherV2GatherV2)Conv1Layer/FGN_0/Tensordot/Shape:output:0(Conv1Layer/FGN_0/Tensordot/free:output:01Conv1Layer/FGN_0/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2%
#Conv1Layer/FGN_0/Tensordot/GatherV2џ
*Conv1Layer/FGN_0/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2,
*Conv1Layer/FGN_0/Tensordot/GatherV2_1/axisг
%Conv1Layer/FGN_0/Tensordot/GatherV2_1GatherV2)Conv1Layer/FGN_0/Tensordot/Shape:output:0(Conv1Layer/FGN_0/Tensordot/axes:output:03Conv1Layer/FGN_0/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2'
%Conv1Layer/FGN_0/Tensordot/GatherV2_1ј
 Conv1Layer/FGN_0/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2"
 Conv1Layer/FGN_0/Tensordot/Const─
Conv1Layer/FGN_0/Tensordot/ProdProd,Conv1Layer/FGN_0/Tensordot/GatherV2:output:0)Conv1Layer/FGN_0/Tensordot/Const:output:0*
T0*
_output_shapes
: 2!
Conv1Layer/FGN_0/Tensordot/Prodњ
"Conv1Layer/FGN_0/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2$
"Conv1Layer/FGN_0/Tensordot/Const_1╠
!Conv1Layer/FGN_0/Tensordot/Prod_1Prod.Conv1Layer/FGN_0/Tensordot/GatherV2_1:output:0+Conv1Layer/FGN_0/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2#
!Conv1Layer/FGN_0/Tensordot/Prod_1њ
&Conv1Layer/FGN_0/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2(
&Conv1Layer/FGN_0/Tensordot/concat/axisЁ
!Conv1Layer/FGN_0/Tensordot/concatConcatV2(Conv1Layer/FGN_0/Tensordot/free:output:0(Conv1Layer/FGN_0/Tensordot/axes:output:0/Conv1Layer/FGN_0/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2#
!Conv1Layer/FGN_0/Tensordot/concatл
 Conv1Layer/FGN_0/Tensordot/stackPack(Conv1Layer/FGN_0/Tensordot/Prod:output:0*Conv1Layer/FGN_0/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2"
 Conv1Layer/FGN_0/Tensordot/stack╔
$Conv1Layer/FGN_0/Tensordot/transpose	Transposeinputs_2*Conv1Layer/FGN_0/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2&
$Conv1Layer/FGN_0/Tensordot/transposeс
"Conv1Layer/FGN_0/Tensordot/ReshapeReshape(Conv1Layer/FGN_0/Tensordot/transpose:y:0)Conv1Layer/FGN_0/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2$
"Conv1Layer/FGN_0/Tensordot/ReshapeР
!Conv1Layer/FGN_0/Tensordot/MatMulMatMul+Conv1Layer/FGN_0/Tensordot/Reshape:output:01Conv1Layer/FGN_0/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2#
!Conv1Layer/FGN_0/Tensordot/MatMulњ
"Conv1Layer/FGN_0/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2$
"Conv1Layer/FGN_0/Tensordot/Const_2ќ
(Conv1Layer/FGN_0/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2*
(Conv1Layer/FGN_0/Tensordot/concat_1/axisњ
#Conv1Layer/FGN_0/Tensordot/concat_1ConcatV2,Conv1Layer/FGN_0/Tensordot/GatherV2:output:0+Conv1Layer/FGN_0/Tensordot/Const_2:output:01Conv1Layer/FGN_0/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2%
#Conv1Layer/FGN_0/Tensordot/concat_1п
Conv1Layer/FGN_0/TensordotReshape+Conv1Layer/FGN_0/Tensordot/MatMul:product:0,Conv1Layer/FGN_0/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_0/Tensordot┐
'Conv1Layer/FGN_0/BiasAdd/ReadVariableOpReadVariableOp0conv1layer_fgn_0_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02)
'Conv1Layer/FGN_0/BiasAdd/ReadVariableOp¤
Conv1Layer/FGN_0/BiasAddBiasAdd#Conv1Layer/FGN_0/Tensordot:output:0/Conv1Layer/FGN_0/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_0/BiasAddЊ
Conv1Layer/FGN_0/ReluRelu!Conv1Layer/FGN_0/BiasAdd:output:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_0/Relu╔
)Conv1Layer/FGN_1/Tensordot/ReadVariableOpReadVariableOp2conv1layer_fgn_1_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02+
)Conv1Layer/FGN_1/Tensordot/ReadVariableOpї
Conv1Layer/FGN_1/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2!
Conv1Layer/FGN_1/Tensordot/axesЌ
Conv1Layer/FGN_1/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv1Layer/FGN_1/Tensordot/freeЌ
 Conv1Layer/FGN_1/Tensordot/ShapeShape#Conv1Layer/FGN_0/Relu:activations:0*
T0*
_output_shapes
:2"
 Conv1Layer/FGN_1/Tensordot/Shapeќ
(Conv1Layer/FGN_1/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2*
(Conv1Layer/FGN_1/Tensordot/GatherV2/axisд
#Conv1Layer/FGN_1/Tensordot/GatherV2GatherV2)Conv1Layer/FGN_1/Tensordot/Shape:output:0(Conv1Layer/FGN_1/Tensordot/free:output:01Conv1Layer/FGN_1/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2%
#Conv1Layer/FGN_1/Tensordot/GatherV2џ
*Conv1Layer/FGN_1/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2,
*Conv1Layer/FGN_1/Tensordot/GatherV2_1/axisг
%Conv1Layer/FGN_1/Tensordot/GatherV2_1GatherV2)Conv1Layer/FGN_1/Tensordot/Shape:output:0(Conv1Layer/FGN_1/Tensordot/axes:output:03Conv1Layer/FGN_1/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2'
%Conv1Layer/FGN_1/Tensordot/GatherV2_1ј
 Conv1Layer/FGN_1/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2"
 Conv1Layer/FGN_1/Tensordot/Const─
Conv1Layer/FGN_1/Tensordot/ProdProd,Conv1Layer/FGN_1/Tensordot/GatherV2:output:0)Conv1Layer/FGN_1/Tensordot/Const:output:0*
T0*
_output_shapes
: 2!
Conv1Layer/FGN_1/Tensordot/Prodњ
"Conv1Layer/FGN_1/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2$
"Conv1Layer/FGN_1/Tensordot/Const_1╠
!Conv1Layer/FGN_1/Tensordot/Prod_1Prod.Conv1Layer/FGN_1/Tensordot/GatherV2_1:output:0+Conv1Layer/FGN_1/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2#
!Conv1Layer/FGN_1/Tensordot/Prod_1њ
&Conv1Layer/FGN_1/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2(
&Conv1Layer/FGN_1/Tensordot/concat/axisЁ
!Conv1Layer/FGN_1/Tensordot/concatConcatV2(Conv1Layer/FGN_1/Tensordot/free:output:0(Conv1Layer/FGN_1/Tensordot/axes:output:0/Conv1Layer/FGN_1/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2#
!Conv1Layer/FGN_1/Tensordot/concatл
 Conv1Layer/FGN_1/Tensordot/stackPack(Conv1Layer/FGN_1/Tensordot/Prod:output:0*Conv1Layer/FGN_1/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2"
 Conv1Layer/FGN_1/Tensordot/stackС
$Conv1Layer/FGN_1/Tensordot/transpose	Transpose#Conv1Layer/FGN_0/Relu:activations:0*Conv1Layer/FGN_1/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2&
$Conv1Layer/FGN_1/Tensordot/transposeс
"Conv1Layer/FGN_1/Tensordot/ReshapeReshape(Conv1Layer/FGN_1/Tensordot/transpose:y:0)Conv1Layer/FGN_1/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2$
"Conv1Layer/FGN_1/Tensordot/ReshapeР
!Conv1Layer/FGN_1/Tensordot/MatMulMatMul+Conv1Layer/FGN_1/Tensordot/Reshape:output:01Conv1Layer/FGN_1/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2#
!Conv1Layer/FGN_1/Tensordot/MatMulњ
"Conv1Layer/FGN_1/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2$
"Conv1Layer/FGN_1/Tensordot/Const_2ќ
(Conv1Layer/FGN_1/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2*
(Conv1Layer/FGN_1/Tensordot/concat_1/axisњ
#Conv1Layer/FGN_1/Tensordot/concat_1ConcatV2,Conv1Layer/FGN_1/Tensordot/GatherV2:output:0+Conv1Layer/FGN_1/Tensordot/Const_2:output:01Conv1Layer/FGN_1/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2%
#Conv1Layer/FGN_1/Tensordot/concat_1п
Conv1Layer/FGN_1/TensordotReshape+Conv1Layer/FGN_1/Tensordot/MatMul:product:0,Conv1Layer/FGN_1/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_1/Tensordot┐
'Conv1Layer/FGN_1/BiasAdd/ReadVariableOpReadVariableOp0conv1layer_fgn_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02)
'Conv1Layer/FGN_1/BiasAdd/ReadVariableOp¤
Conv1Layer/FGN_1/BiasAddBiasAdd#Conv1Layer/FGN_1/Tensordot:output:0/Conv1Layer/FGN_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_1/BiasAddЊ
Conv1Layer/FGN_1/ReluRelu!Conv1Layer/FGN_1/BiasAdd:output:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_1/Relu¤
+Conv1Layer/FGN_out/Tensordot/ReadVariableOpReadVariableOp4conv1layer_fgn_out_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02-
+Conv1Layer/FGN_out/Tensordot/ReadVariableOpљ
!Conv1Layer/FGN_out/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2#
!Conv1Layer/FGN_out/Tensordot/axesЏ
!Conv1Layer/FGN_out/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2#
!Conv1Layer/FGN_out/Tensordot/freeЏ
"Conv1Layer/FGN_out/Tensordot/ShapeShape#Conv1Layer/FGN_1/Relu:activations:0*
T0*
_output_shapes
:2$
"Conv1Layer/FGN_out/Tensordot/Shapeџ
*Conv1Layer/FGN_out/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2,
*Conv1Layer/FGN_out/Tensordot/GatherV2/axis░
%Conv1Layer/FGN_out/Tensordot/GatherV2GatherV2+Conv1Layer/FGN_out/Tensordot/Shape:output:0*Conv1Layer/FGN_out/Tensordot/free:output:03Conv1Layer/FGN_out/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2'
%Conv1Layer/FGN_out/Tensordot/GatherV2ъ
,Conv1Layer/FGN_out/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2.
,Conv1Layer/FGN_out/Tensordot/GatherV2_1/axisХ
'Conv1Layer/FGN_out/Tensordot/GatherV2_1GatherV2+Conv1Layer/FGN_out/Tensordot/Shape:output:0*Conv1Layer/FGN_out/Tensordot/axes:output:05Conv1Layer/FGN_out/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2)
'Conv1Layer/FGN_out/Tensordot/GatherV2_1њ
"Conv1Layer/FGN_out/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2$
"Conv1Layer/FGN_out/Tensordot/Const╠
!Conv1Layer/FGN_out/Tensordot/ProdProd.Conv1Layer/FGN_out/Tensordot/GatherV2:output:0+Conv1Layer/FGN_out/Tensordot/Const:output:0*
T0*
_output_shapes
: 2#
!Conv1Layer/FGN_out/Tensordot/Prodќ
$Conv1Layer/FGN_out/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2&
$Conv1Layer/FGN_out/Tensordot/Const_1н
#Conv1Layer/FGN_out/Tensordot/Prod_1Prod0Conv1Layer/FGN_out/Tensordot/GatherV2_1:output:0-Conv1Layer/FGN_out/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2%
#Conv1Layer/FGN_out/Tensordot/Prod_1ќ
(Conv1Layer/FGN_out/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2*
(Conv1Layer/FGN_out/Tensordot/concat/axisЈ
#Conv1Layer/FGN_out/Tensordot/concatConcatV2*Conv1Layer/FGN_out/Tensordot/free:output:0*Conv1Layer/FGN_out/Tensordot/axes:output:01Conv1Layer/FGN_out/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2%
#Conv1Layer/FGN_out/Tensordot/concatп
"Conv1Layer/FGN_out/Tensordot/stackPack*Conv1Layer/FGN_out/Tensordot/Prod:output:0,Conv1Layer/FGN_out/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2$
"Conv1Layer/FGN_out/Tensordot/stackЖ
&Conv1Layer/FGN_out/Tensordot/transpose	Transpose#Conv1Layer/FGN_1/Relu:activations:0,Conv1Layer/FGN_out/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2(
&Conv1Layer/FGN_out/Tensordot/transposeв
$Conv1Layer/FGN_out/Tensordot/ReshapeReshape*Conv1Layer/FGN_out/Tensordot/transpose:y:0+Conv1Layer/FGN_out/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2&
$Conv1Layer/FGN_out/Tensordot/ReshapeЖ
#Conv1Layer/FGN_out/Tensordot/MatMulMatMul-Conv1Layer/FGN_out/Tensordot/Reshape:output:03Conv1Layer/FGN_out/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2%
#Conv1Layer/FGN_out/Tensordot/MatMulќ
$Conv1Layer/FGN_out/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2&
$Conv1Layer/FGN_out/Tensordot/Const_2џ
*Conv1Layer/FGN_out/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2,
*Conv1Layer/FGN_out/Tensordot/concat_1/axisю
%Conv1Layer/FGN_out/Tensordot/concat_1ConcatV2.Conv1Layer/FGN_out/Tensordot/GatherV2:output:0-Conv1Layer/FGN_out/Tensordot/Const_2:output:03Conv1Layer/FGN_out/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2'
%Conv1Layer/FGN_out/Tensordot/concat_1Я
Conv1Layer/FGN_out/TensordotReshape-Conv1Layer/FGN_out/Tensordot/MatMul:product:0.Conv1Layer/FGN_out/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_out/Tensordot┼
)Conv1Layer/FGN_out/BiasAdd/ReadVariableOpReadVariableOp2conv1layer_fgn_out_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)Conv1Layer/FGN_out/BiasAdd/ReadVariableOpО
Conv1Layer/FGN_out/BiasAddBiasAdd%Conv1Layer/FGN_out/Tensordot:output:01Conv1Layer/FGN_out/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_out/BiasAddЃ
Conv1Layer/Reshape/shape/0Const*
_output_shapes
: *
dtype0*
valueB :
         2
Conv1Layer/Reshape/shape/0z
Conv1Layer/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :2
Conv1Layer/Reshape/shape/3z
Conv1Layer/Reshape/shape/4Const*
_output_shapes
: *
dtype0*
value	B :2
Conv1Layer/Reshape/shape/4Ъ
Conv1Layer/Reshape/shapePack#Conv1Layer/Reshape/shape/0:output:0!Conv1Layer/strided_slice:output:0!Conv1Layer/strided_slice:output:0#Conv1Layer/Reshape/shape/3:output:0#Conv1Layer/Reshape/shape/4:output:0*
N*
T0*
_output_shapes
:2
Conv1Layer/Reshape/shape╦
Conv1Layer/ReshapeReshape#Conv1Layer/FGN_out/BiasAdd:output:0!Conv1Layer/Reshape/shape:output:0*
T0*E
_output_shapes3
1:/                           2
Conv1Layer/ReshapeЎ
 Conv1Layer/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2"
 Conv1Layer/strided_slice_1/stackЮ
"Conv1Layer/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"            2$
"Conv1Layer/strided_slice_1/stack_1Ю
"Conv1Layer/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2$
"Conv1Layer/strided_slice_1/stack_2╠
Conv1Layer/strided_slice_1StridedSliceinputs_1)Conv1Layer/strided_slice_1/stack:output:0+Conv1Layer/strided_slice_1/stack_1:output:0+Conv1Layer/strided_slice_1/stack_2:output:0*
Index0*
T0*3
_output_shapes!
:         *
ellipsis_mask*
new_axis_mask2
Conv1Layer/strided_slice_1Д
Conv1Layer/mulMulConv1Layer/Reshape:output:0#Conv1Layer/strided_slice_1:output:0*
T0*3
_output_shapes!
:         2
Conv1Layer/mul╣
Conv1Layer/einsum/EinsumEinsumConv1Layer/mul:z:0inputs_0*
N*
T0*+
_output_shapes
:         *
equationabicf,aif->abc2
Conv1Layer/einsum/Einsum«
 Conv1Layer/MatMul/ReadVariableOpReadVariableOp)conv1layer_matmul_readvariableop_resource*
_output_shapes

:*
dtype02"
 Conv1Layer/MatMul/ReadVariableOpА
Conv1Layer/MatMulBatchMatMulV2inputs_0(Conv1Layer/MatMul/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
Conv1Layer/MatMulъ
Conv1Layer/addAddV2!Conv1Layer/einsum/Einsum:output:0Conv1Layer/MatMul:output:0*
T0*+
_output_shapes
:         2
Conv1Layer/addГ
!Conv1Layer/BiasAdd/ReadVariableOpReadVariableOp*conv1layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02#
!Conv1Layer/BiasAdd/ReadVariableOpе
Conv1Layer/BiasAddBiasAddConv1Layer/add:z:0)Conv1Layer/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
Conv1Layer/BiasAddz
Conv1Layer/EluEluConv1Layer/BiasAdd:output:0*
T0*+
_output_shapes
:         2
Conv1Layer/Elu║
$OutputLayer/Tensordot/ReadVariableOpReadVariableOp-outputlayer_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02&
$OutputLayer/Tensordot/ReadVariableOpѓ
OutputLayer/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2
OutputLayer/Tensordot/axesЅ
OutputLayer/Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB"       2
OutputLayer/Tensordot/freeє
OutputLayer/Tensordot/ShapeShapeConv1Layer/Elu:activations:0*
T0*
_output_shapes
:2
OutputLayer/Tensordot/Shapeї
#OutputLayer/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2%
#OutputLayer/Tensordot/GatherV2/axisЇ
OutputLayer/Tensordot/GatherV2GatherV2$OutputLayer/Tensordot/Shape:output:0#OutputLayer/Tensordot/free:output:0,OutputLayer/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2 
OutputLayer/Tensordot/GatherV2љ
%OutputLayer/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2'
%OutputLayer/Tensordot/GatherV2_1/axisЊ
 OutputLayer/Tensordot/GatherV2_1GatherV2$OutputLayer/Tensordot/Shape:output:0#OutputLayer/Tensordot/axes:output:0.OutputLayer/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2"
 OutputLayer/Tensordot/GatherV2_1ё
OutputLayer/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
OutputLayer/Tensordot/Const░
OutputLayer/Tensordot/ProdProd'OutputLayer/Tensordot/GatherV2:output:0$OutputLayer/Tensordot/Const:output:0*
T0*
_output_shapes
: 2
OutputLayer/Tensordot/Prodѕ
OutputLayer/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2
OutputLayer/Tensordot/Const_1И
OutputLayer/Tensordot/Prod_1Prod)OutputLayer/Tensordot/GatherV2_1:output:0&OutputLayer/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2
OutputLayer/Tensordot/Prod_1ѕ
!OutputLayer/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2#
!OutputLayer/Tensordot/concat/axisВ
OutputLayer/Tensordot/concatConcatV2#OutputLayer/Tensordot/free:output:0#OutputLayer/Tensordot/axes:output:0*OutputLayer/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2
OutputLayer/Tensordot/concat╝
OutputLayer/Tensordot/stackPack#OutputLayer/Tensordot/Prod:output:0%OutputLayer/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2
OutputLayer/Tensordot/stack╩
OutputLayer/Tensordot/transpose	TransposeConv1Layer/Elu:activations:0%OutputLayer/Tensordot/concat:output:0*
T0*+
_output_shapes
:         2!
OutputLayer/Tensordot/transpose¤
OutputLayer/Tensordot/ReshapeReshape#OutputLayer/Tensordot/transpose:y:0$OutputLayer/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2
OutputLayer/Tensordot/Reshape╬
OutputLayer/Tensordot/MatMulMatMul&OutputLayer/Tensordot/Reshape:output:0,OutputLayer/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
OutputLayer/Tensordot/MatMulѕ
OutputLayer/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2
OutputLayer/Tensordot/Const_2ї
#OutputLayer/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2%
#OutputLayer/Tensordot/concat_1/axisщ
OutputLayer/Tensordot/concat_1ConcatV2'OutputLayer/Tensordot/GatherV2:output:0&OutputLayer/Tensordot/Const_2:output:0,OutputLayer/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2 
OutputLayer/Tensordot/concat_1└
OutputLayer/TensordotReshape&OutputLayer/Tensordot/MatMul:product:0'OutputLayer/Tensordot/concat_1:output:0*
T0*+
_output_shapes
:         2
OutputLayer/Tensordot░
"OutputLayer/BiasAdd/ReadVariableOpReadVariableOp+outputlayer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02$
"OutputLayer/BiasAdd/ReadVariableOpи
OutputLayer/BiasAddBiasAddOutputLayer/Tensordot:output:0*OutputLayer/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
OutputLayer/BiasAdd}
OutputLayer/EluEluOutputLayer/BiasAdd:output:0*
T0*+
_output_shapes
:         2
OutputLayer/Eluu
IdentityIdentityOutputLayer/Elu:activations:0*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         :::::::::::U Q
+
_output_shapes
:         
"
_user_specified_name
inputs/0:UQ
+
_output_shapes
:         
"
_user_specified_name
inputs/1:YU
/
_output_shapes
:         
"
_user_specified_name
inputs/2:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
Ў└
Ј
A__inference_MyModel_layer_call_and_return_conditional_losses_1028
inputs_0
inputs_1
inputs_26
2conv1layer_fgn_0_tensordot_readvariableop_resource4
0conv1layer_fgn_0_biasadd_readvariableop_resource6
2conv1layer_fgn_1_tensordot_readvariableop_resource4
0conv1layer_fgn_1_biasadd_readvariableop_resource8
4conv1layer_fgn_out_tensordot_readvariableop_resource6
2conv1layer_fgn_out_biasadd_readvariableop_resource-
)conv1layer_matmul_readvariableop_resource.
*conv1layer_biasadd_readvariableop_resource1
-outputlayer_tensordot_readvariableop_resource/
+outputlayer_biasadd_readvariableop_resource
identityѕ\
Conv1Layer/ShapeShapeinputs_0*
T0*
_output_shapes
:2
Conv1Layer/ShapeЊ
Conv1Layer/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2 
Conv1Layer/strided_slice/stackЌ
 Conv1Layer/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2"
 Conv1Layer/strided_slice/stack_1ј
 Conv1Layer/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2"
 Conv1Layer/strided_slice/stack_2ц
Conv1Layer/strided_sliceStridedSliceConv1Layer/Shape:output:0'Conv1Layer/strided_slice/stack:output:0)Conv1Layer/strided_slice/stack_1:output:0)Conv1Layer/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
Conv1Layer/strided_slice╔
)Conv1Layer/FGN_0/Tensordot/ReadVariableOpReadVariableOp2conv1layer_fgn_0_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02+
)Conv1Layer/FGN_0/Tensordot/ReadVariableOpї
Conv1Layer/FGN_0/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2!
Conv1Layer/FGN_0/Tensordot/axesЌ
Conv1Layer/FGN_0/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv1Layer/FGN_0/Tensordot/free|
 Conv1Layer/FGN_0/Tensordot/ShapeShapeinputs_2*
T0*
_output_shapes
:2"
 Conv1Layer/FGN_0/Tensordot/Shapeќ
(Conv1Layer/FGN_0/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2*
(Conv1Layer/FGN_0/Tensordot/GatherV2/axisд
#Conv1Layer/FGN_0/Tensordot/GatherV2GatherV2)Conv1Layer/FGN_0/Tensordot/Shape:output:0(Conv1Layer/FGN_0/Tensordot/free:output:01Conv1Layer/FGN_0/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2%
#Conv1Layer/FGN_0/Tensordot/GatherV2џ
*Conv1Layer/FGN_0/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2,
*Conv1Layer/FGN_0/Tensordot/GatherV2_1/axisг
%Conv1Layer/FGN_0/Tensordot/GatherV2_1GatherV2)Conv1Layer/FGN_0/Tensordot/Shape:output:0(Conv1Layer/FGN_0/Tensordot/axes:output:03Conv1Layer/FGN_0/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2'
%Conv1Layer/FGN_0/Tensordot/GatherV2_1ј
 Conv1Layer/FGN_0/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2"
 Conv1Layer/FGN_0/Tensordot/Const─
Conv1Layer/FGN_0/Tensordot/ProdProd,Conv1Layer/FGN_0/Tensordot/GatherV2:output:0)Conv1Layer/FGN_0/Tensordot/Const:output:0*
T0*
_output_shapes
: 2!
Conv1Layer/FGN_0/Tensordot/Prodњ
"Conv1Layer/FGN_0/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2$
"Conv1Layer/FGN_0/Tensordot/Const_1╠
!Conv1Layer/FGN_0/Tensordot/Prod_1Prod.Conv1Layer/FGN_0/Tensordot/GatherV2_1:output:0+Conv1Layer/FGN_0/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2#
!Conv1Layer/FGN_0/Tensordot/Prod_1њ
&Conv1Layer/FGN_0/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2(
&Conv1Layer/FGN_0/Tensordot/concat/axisЁ
!Conv1Layer/FGN_0/Tensordot/concatConcatV2(Conv1Layer/FGN_0/Tensordot/free:output:0(Conv1Layer/FGN_0/Tensordot/axes:output:0/Conv1Layer/FGN_0/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2#
!Conv1Layer/FGN_0/Tensordot/concatл
 Conv1Layer/FGN_0/Tensordot/stackPack(Conv1Layer/FGN_0/Tensordot/Prod:output:0*Conv1Layer/FGN_0/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2"
 Conv1Layer/FGN_0/Tensordot/stack╔
$Conv1Layer/FGN_0/Tensordot/transpose	Transposeinputs_2*Conv1Layer/FGN_0/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2&
$Conv1Layer/FGN_0/Tensordot/transposeс
"Conv1Layer/FGN_0/Tensordot/ReshapeReshape(Conv1Layer/FGN_0/Tensordot/transpose:y:0)Conv1Layer/FGN_0/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2$
"Conv1Layer/FGN_0/Tensordot/ReshapeР
!Conv1Layer/FGN_0/Tensordot/MatMulMatMul+Conv1Layer/FGN_0/Tensordot/Reshape:output:01Conv1Layer/FGN_0/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2#
!Conv1Layer/FGN_0/Tensordot/MatMulњ
"Conv1Layer/FGN_0/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2$
"Conv1Layer/FGN_0/Tensordot/Const_2ќ
(Conv1Layer/FGN_0/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2*
(Conv1Layer/FGN_0/Tensordot/concat_1/axisњ
#Conv1Layer/FGN_0/Tensordot/concat_1ConcatV2,Conv1Layer/FGN_0/Tensordot/GatherV2:output:0+Conv1Layer/FGN_0/Tensordot/Const_2:output:01Conv1Layer/FGN_0/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2%
#Conv1Layer/FGN_0/Tensordot/concat_1п
Conv1Layer/FGN_0/TensordotReshape+Conv1Layer/FGN_0/Tensordot/MatMul:product:0,Conv1Layer/FGN_0/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_0/Tensordot┐
'Conv1Layer/FGN_0/BiasAdd/ReadVariableOpReadVariableOp0conv1layer_fgn_0_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02)
'Conv1Layer/FGN_0/BiasAdd/ReadVariableOp¤
Conv1Layer/FGN_0/BiasAddBiasAdd#Conv1Layer/FGN_0/Tensordot:output:0/Conv1Layer/FGN_0/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_0/BiasAddЊ
Conv1Layer/FGN_0/ReluRelu!Conv1Layer/FGN_0/BiasAdd:output:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_0/Relu╔
)Conv1Layer/FGN_1/Tensordot/ReadVariableOpReadVariableOp2conv1layer_fgn_1_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02+
)Conv1Layer/FGN_1/Tensordot/ReadVariableOpї
Conv1Layer/FGN_1/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2!
Conv1Layer/FGN_1/Tensordot/axesЌ
Conv1Layer/FGN_1/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2!
Conv1Layer/FGN_1/Tensordot/freeЌ
 Conv1Layer/FGN_1/Tensordot/ShapeShape#Conv1Layer/FGN_0/Relu:activations:0*
T0*
_output_shapes
:2"
 Conv1Layer/FGN_1/Tensordot/Shapeќ
(Conv1Layer/FGN_1/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2*
(Conv1Layer/FGN_1/Tensordot/GatherV2/axisд
#Conv1Layer/FGN_1/Tensordot/GatherV2GatherV2)Conv1Layer/FGN_1/Tensordot/Shape:output:0(Conv1Layer/FGN_1/Tensordot/free:output:01Conv1Layer/FGN_1/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2%
#Conv1Layer/FGN_1/Tensordot/GatherV2џ
*Conv1Layer/FGN_1/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2,
*Conv1Layer/FGN_1/Tensordot/GatherV2_1/axisг
%Conv1Layer/FGN_1/Tensordot/GatherV2_1GatherV2)Conv1Layer/FGN_1/Tensordot/Shape:output:0(Conv1Layer/FGN_1/Tensordot/axes:output:03Conv1Layer/FGN_1/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2'
%Conv1Layer/FGN_1/Tensordot/GatherV2_1ј
 Conv1Layer/FGN_1/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2"
 Conv1Layer/FGN_1/Tensordot/Const─
Conv1Layer/FGN_1/Tensordot/ProdProd,Conv1Layer/FGN_1/Tensordot/GatherV2:output:0)Conv1Layer/FGN_1/Tensordot/Const:output:0*
T0*
_output_shapes
: 2!
Conv1Layer/FGN_1/Tensordot/Prodњ
"Conv1Layer/FGN_1/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2$
"Conv1Layer/FGN_1/Tensordot/Const_1╠
!Conv1Layer/FGN_1/Tensordot/Prod_1Prod.Conv1Layer/FGN_1/Tensordot/GatherV2_1:output:0+Conv1Layer/FGN_1/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2#
!Conv1Layer/FGN_1/Tensordot/Prod_1њ
&Conv1Layer/FGN_1/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2(
&Conv1Layer/FGN_1/Tensordot/concat/axisЁ
!Conv1Layer/FGN_1/Tensordot/concatConcatV2(Conv1Layer/FGN_1/Tensordot/free:output:0(Conv1Layer/FGN_1/Tensordot/axes:output:0/Conv1Layer/FGN_1/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2#
!Conv1Layer/FGN_1/Tensordot/concatл
 Conv1Layer/FGN_1/Tensordot/stackPack(Conv1Layer/FGN_1/Tensordot/Prod:output:0*Conv1Layer/FGN_1/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2"
 Conv1Layer/FGN_1/Tensordot/stackС
$Conv1Layer/FGN_1/Tensordot/transpose	Transpose#Conv1Layer/FGN_0/Relu:activations:0*Conv1Layer/FGN_1/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2&
$Conv1Layer/FGN_1/Tensordot/transposeс
"Conv1Layer/FGN_1/Tensordot/ReshapeReshape(Conv1Layer/FGN_1/Tensordot/transpose:y:0)Conv1Layer/FGN_1/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2$
"Conv1Layer/FGN_1/Tensordot/ReshapeР
!Conv1Layer/FGN_1/Tensordot/MatMulMatMul+Conv1Layer/FGN_1/Tensordot/Reshape:output:01Conv1Layer/FGN_1/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2#
!Conv1Layer/FGN_1/Tensordot/MatMulњ
"Conv1Layer/FGN_1/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2$
"Conv1Layer/FGN_1/Tensordot/Const_2ќ
(Conv1Layer/FGN_1/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2*
(Conv1Layer/FGN_1/Tensordot/concat_1/axisњ
#Conv1Layer/FGN_1/Tensordot/concat_1ConcatV2,Conv1Layer/FGN_1/Tensordot/GatherV2:output:0+Conv1Layer/FGN_1/Tensordot/Const_2:output:01Conv1Layer/FGN_1/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2%
#Conv1Layer/FGN_1/Tensordot/concat_1п
Conv1Layer/FGN_1/TensordotReshape+Conv1Layer/FGN_1/Tensordot/MatMul:product:0,Conv1Layer/FGN_1/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_1/Tensordot┐
'Conv1Layer/FGN_1/BiasAdd/ReadVariableOpReadVariableOp0conv1layer_fgn_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02)
'Conv1Layer/FGN_1/BiasAdd/ReadVariableOp¤
Conv1Layer/FGN_1/BiasAddBiasAdd#Conv1Layer/FGN_1/Tensordot:output:0/Conv1Layer/FGN_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_1/BiasAddЊ
Conv1Layer/FGN_1/ReluRelu!Conv1Layer/FGN_1/BiasAdd:output:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_1/Relu¤
+Conv1Layer/FGN_out/Tensordot/ReadVariableOpReadVariableOp4conv1layer_fgn_out_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02-
+Conv1Layer/FGN_out/Tensordot/ReadVariableOpљ
!Conv1Layer/FGN_out/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2#
!Conv1Layer/FGN_out/Tensordot/axesЏ
!Conv1Layer/FGN_out/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2#
!Conv1Layer/FGN_out/Tensordot/freeЏ
"Conv1Layer/FGN_out/Tensordot/ShapeShape#Conv1Layer/FGN_1/Relu:activations:0*
T0*
_output_shapes
:2$
"Conv1Layer/FGN_out/Tensordot/Shapeџ
*Conv1Layer/FGN_out/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2,
*Conv1Layer/FGN_out/Tensordot/GatherV2/axis░
%Conv1Layer/FGN_out/Tensordot/GatherV2GatherV2+Conv1Layer/FGN_out/Tensordot/Shape:output:0*Conv1Layer/FGN_out/Tensordot/free:output:03Conv1Layer/FGN_out/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2'
%Conv1Layer/FGN_out/Tensordot/GatherV2ъ
,Conv1Layer/FGN_out/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2.
,Conv1Layer/FGN_out/Tensordot/GatherV2_1/axisХ
'Conv1Layer/FGN_out/Tensordot/GatherV2_1GatherV2+Conv1Layer/FGN_out/Tensordot/Shape:output:0*Conv1Layer/FGN_out/Tensordot/axes:output:05Conv1Layer/FGN_out/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2)
'Conv1Layer/FGN_out/Tensordot/GatherV2_1њ
"Conv1Layer/FGN_out/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2$
"Conv1Layer/FGN_out/Tensordot/Const╠
!Conv1Layer/FGN_out/Tensordot/ProdProd.Conv1Layer/FGN_out/Tensordot/GatherV2:output:0+Conv1Layer/FGN_out/Tensordot/Const:output:0*
T0*
_output_shapes
: 2#
!Conv1Layer/FGN_out/Tensordot/Prodќ
$Conv1Layer/FGN_out/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2&
$Conv1Layer/FGN_out/Tensordot/Const_1н
#Conv1Layer/FGN_out/Tensordot/Prod_1Prod0Conv1Layer/FGN_out/Tensordot/GatherV2_1:output:0-Conv1Layer/FGN_out/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2%
#Conv1Layer/FGN_out/Tensordot/Prod_1ќ
(Conv1Layer/FGN_out/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2*
(Conv1Layer/FGN_out/Tensordot/concat/axisЈ
#Conv1Layer/FGN_out/Tensordot/concatConcatV2*Conv1Layer/FGN_out/Tensordot/free:output:0*Conv1Layer/FGN_out/Tensordot/axes:output:01Conv1Layer/FGN_out/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2%
#Conv1Layer/FGN_out/Tensordot/concatп
"Conv1Layer/FGN_out/Tensordot/stackPack*Conv1Layer/FGN_out/Tensordot/Prod:output:0,Conv1Layer/FGN_out/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2$
"Conv1Layer/FGN_out/Tensordot/stackЖ
&Conv1Layer/FGN_out/Tensordot/transpose	Transpose#Conv1Layer/FGN_1/Relu:activations:0,Conv1Layer/FGN_out/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2(
&Conv1Layer/FGN_out/Tensordot/transposeв
$Conv1Layer/FGN_out/Tensordot/ReshapeReshape*Conv1Layer/FGN_out/Tensordot/transpose:y:0+Conv1Layer/FGN_out/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2&
$Conv1Layer/FGN_out/Tensordot/ReshapeЖ
#Conv1Layer/FGN_out/Tensordot/MatMulMatMul-Conv1Layer/FGN_out/Tensordot/Reshape:output:03Conv1Layer/FGN_out/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2%
#Conv1Layer/FGN_out/Tensordot/MatMulќ
$Conv1Layer/FGN_out/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2&
$Conv1Layer/FGN_out/Tensordot/Const_2џ
*Conv1Layer/FGN_out/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2,
*Conv1Layer/FGN_out/Tensordot/concat_1/axisю
%Conv1Layer/FGN_out/Tensordot/concat_1ConcatV2.Conv1Layer/FGN_out/Tensordot/GatherV2:output:0-Conv1Layer/FGN_out/Tensordot/Const_2:output:03Conv1Layer/FGN_out/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2'
%Conv1Layer/FGN_out/Tensordot/concat_1Я
Conv1Layer/FGN_out/TensordotReshape-Conv1Layer/FGN_out/Tensordot/MatMul:product:0.Conv1Layer/FGN_out/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_out/Tensordot┼
)Conv1Layer/FGN_out/BiasAdd/ReadVariableOpReadVariableOp2conv1layer_fgn_out_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)Conv1Layer/FGN_out/BiasAdd/ReadVariableOpО
Conv1Layer/FGN_out/BiasAddBiasAdd%Conv1Layer/FGN_out/Tensordot:output:01Conv1Layer/FGN_out/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
Conv1Layer/FGN_out/BiasAddЃ
Conv1Layer/Reshape/shape/0Const*
_output_shapes
: *
dtype0*
valueB :
         2
Conv1Layer/Reshape/shape/0z
Conv1Layer/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :2
Conv1Layer/Reshape/shape/3z
Conv1Layer/Reshape/shape/4Const*
_output_shapes
: *
dtype0*
value	B :2
Conv1Layer/Reshape/shape/4Ъ
Conv1Layer/Reshape/shapePack#Conv1Layer/Reshape/shape/0:output:0!Conv1Layer/strided_slice:output:0!Conv1Layer/strided_slice:output:0#Conv1Layer/Reshape/shape/3:output:0#Conv1Layer/Reshape/shape/4:output:0*
N*
T0*
_output_shapes
:2
Conv1Layer/Reshape/shape╦
Conv1Layer/ReshapeReshape#Conv1Layer/FGN_out/BiasAdd:output:0!Conv1Layer/Reshape/shape:output:0*
T0*E
_output_shapes3
1:/                           2
Conv1Layer/ReshapeЎ
 Conv1Layer/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2"
 Conv1Layer/strided_slice_1/stackЮ
"Conv1Layer/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"            2$
"Conv1Layer/strided_slice_1/stack_1Ю
"Conv1Layer/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2$
"Conv1Layer/strided_slice_1/stack_2╠
Conv1Layer/strided_slice_1StridedSliceinputs_1)Conv1Layer/strided_slice_1/stack:output:0+Conv1Layer/strided_slice_1/stack_1:output:0+Conv1Layer/strided_slice_1/stack_2:output:0*
Index0*
T0*3
_output_shapes!
:         *
ellipsis_mask*
new_axis_mask2
Conv1Layer/strided_slice_1Д
Conv1Layer/mulMulConv1Layer/Reshape:output:0#Conv1Layer/strided_slice_1:output:0*
T0*3
_output_shapes!
:         2
Conv1Layer/mul╣
Conv1Layer/einsum/EinsumEinsumConv1Layer/mul:z:0inputs_0*
N*
T0*+
_output_shapes
:         *
equationabicf,aif->abc2
Conv1Layer/einsum/Einsum«
 Conv1Layer/MatMul/ReadVariableOpReadVariableOp)conv1layer_matmul_readvariableop_resource*
_output_shapes

:*
dtype02"
 Conv1Layer/MatMul/ReadVariableOpА
Conv1Layer/MatMulBatchMatMulV2inputs_0(Conv1Layer/MatMul/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
Conv1Layer/MatMulъ
Conv1Layer/addAddV2!Conv1Layer/einsum/Einsum:output:0Conv1Layer/MatMul:output:0*
T0*+
_output_shapes
:         2
Conv1Layer/addГ
!Conv1Layer/BiasAdd/ReadVariableOpReadVariableOp*conv1layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02#
!Conv1Layer/BiasAdd/ReadVariableOpе
Conv1Layer/BiasAddBiasAddConv1Layer/add:z:0)Conv1Layer/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
Conv1Layer/BiasAddz
Conv1Layer/EluEluConv1Layer/BiasAdd:output:0*
T0*+
_output_shapes
:         2
Conv1Layer/Elu║
$OutputLayer/Tensordot/ReadVariableOpReadVariableOp-outputlayer_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02&
$OutputLayer/Tensordot/ReadVariableOpѓ
OutputLayer/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2
OutputLayer/Tensordot/axesЅ
OutputLayer/Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB"       2
OutputLayer/Tensordot/freeє
OutputLayer/Tensordot/ShapeShapeConv1Layer/Elu:activations:0*
T0*
_output_shapes
:2
OutputLayer/Tensordot/Shapeї
#OutputLayer/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2%
#OutputLayer/Tensordot/GatherV2/axisЇ
OutputLayer/Tensordot/GatherV2GatherV2$OutputLayer/Tensordot/Shape:output:0#OutputLayer/Tensordot/free:output:0,OutputLayer/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2 
OutputLayer/Tensordot/GatherV2љ
%OutputLayer/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2'
%OutputLayer/Tensordot/GatherV2_1/axisЊ
 OutputLayer/Tensordot/GatherV2_1GatherV2$OutputLayer/Tensordot/Shape:output:0#OutputLayer/Tensordot/axes:output:0.OutputLayer/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2"
 OutputLayer/Tensordot/GatherV2_1ё
OutputLayer/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
OutputLayer/Tensordot/Const░
OutputLayer/Tensordot/ProdProd'OutputLayer/Tensordot/GatherV2:output:0$OutputLayer/Tensordot/Const:output:0*
T0*
_output_shapes
: 2
OutputLayer/Tensordot/Prodѕ
OutputLayer/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2
OutputLayer/Tensordot/Const_1И
OutputLayer/Tensordot/Prod_1Prod)OutputLayer/Tensordot/GatherV2_1:output:0&OutputLayer/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2
OutputLayer/Tensordot/Prod_1ѕ
!OutputLayer/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2#
!OutputLayer/Tensordot/concat/axisВ
OutputLayer/Tensordot/concatConcatV2#OutputLayer/Tensordot/free:output:0#OutputLayer/Tensordot/axes:output:0*OutputLayer/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2
OutputLayer/Tensordot/concat╝
OutputLayer/Tensordot/stackPack#OutputLayer/Tensordot/Prod:output:0%OutputLayer/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2
OutputLayer/Tensordot/stack╩
OutputLayer/Tensordot/transpose	TransposeConv1Layer/Elu:activations:0%OutputLayer/Tensordot/concat:output:0*
T0*+
_output_shapes
:         2!
OutputLayer/Tensordot/transpose¤
OutputLayer/Tensordot/ReshapeReshape#OutputLayer/Tensordot/transpose:y:0$OutputLayer/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2
OutputLayer/Tensordot/Reshape╬
OutputLayer/Tensordot/MatMulMatMul&OutputLayer/Tensordot/Reshape:output:0,OutputLayer/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
OutputLayer/Tensordot/MatMulѕ
OutputLayer/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2
OutputLayer/Tensordot/Const_2ї
#OutputLayer/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2%
#OutputLayer/Tensordot/concat_1/axisщ
OutputLayer/Tensordot/concat_1ConcatV2'OutputLayer/Tensordot/GatherV2:output:0&OutputLayer/Tensordot/Const_2:output:0,OutputLayer/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2 
OutputLayer/Tensordot/concat_1└
OutputLayer/TensordotReshape&OutputLayer/Tensordot/MatMul:product:0'OutputLayer/Tensordot/concat_1:output:0*
T0*+
_output_shapes
:         2
OutputLayer/Tensordot░
"OutputLayer/BiasAdd/ReadVariableOpReadVariableOp+outputlayer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02$
"OutputLayer/BiasAdd/ReadVariableOpи
OutputLayer/BiasAddBiasAddOutputLayer/Tensordot:output:0*OutputLayer/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
OutputLayer/BiasAdd}
OutputLayer/EluEluOutputLayer/BiasAdd:output:0*
T0*+
_output_shapes
:         2
OutputLayer/Eluu
IdentityIdentityOutputLayer/Elu:activations:0*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         :::::::::::U Q
+
_output_shapes
:         
"
_user_specified_name
inputs/0:UQ
+
_output_shapes
:         
"
_user_specified_name
inputs/1:YU
/
_output_shapes
:         
"
_user_specified_name
inputs/2:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
џђ
о
D__inference_Conv1Layer_layer_call_and_return_conditional_losses_1329
inputs_0
inputs_1
inputs_2+
'fgn_0_tensordot_readvariableop_resource)
%fgn_0_biasadd_readvariableop_resource+
'fgn_1_tensordot_readvariableop_resource)
%fgn_1_biasadd_readvariableop_resource-
)fgn_out_tensordot_readvariableop_resource+
'fgn_out_biasadd_readvariableop_resource"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityѕF
ShapeShapeinputs_0*
T0*
_output_shapes
:2
Shape}
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2
strided_slice/stackЂ
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2Р
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_sliceе
FGN_0/Tensordot/ReadVariableOpReadVariableOp'fgn_0_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02 
FGN_0/Tensordot/ReadVariableOpv
FGN_0/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2
FGN_0/Tensordot/axesЂ
FGN_0/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2
FGN_0/Tensordot/freef
FGN_0/Tensordot/ShapeShapeinputs_2*
T0*
_output_shapes
:2
FGN_0/Tensordot/Shapeђ
FGN_0/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_0/Tensordot/GatherV2/axis№
FGN_0/Tensordot/GatherV2GatherV2FGN_0/Tensordot/Shape:output:0FGN_0/Tensordot/free:output:0&FGN_0/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_0/Tensordot/GatherV2ё
FGN_0/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2!
FGN_0/Tensordot/GatherV2_1/axisш
FGN_0/Tensordot/GatherV2_1GatherV2FGN_0/Tensordot/Shape:output:0FGN_0/Tensordot/axes:output:0(FGN_0/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_0/Tensordot/GatherV2_1x
FGN_0/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
FGN_0/Tensordot/Constў
FGN_0/Tensordot/ProdProd!FGN_0/Tensordot/GatherV2:output:0FGN_0/Tensordot/Const:output:0*
T0*
_output_shapes
: 2
FGN_0/Tensordot/Prod|
FGN_0/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2
FGN_0/Tensordot/Const_1а
FGN_0/Tensordot/Prod_1Prod#FGN_0/Tensordot/GatherV2_1:output:0 FGN_0/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2
FGN_0/Tensordot/Prod_1|
FGN_0/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_0/Tensordot/concat/axis╬
FGN_0/Tensordot/concatConcatV2FGN_0/Tensordot/free:output:0FGN_0/Tensordot/axes:output:0$FGN_0/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_0/Tensordot/concatц
FGN_0/Tensordot/stackPackFGN_0/Tensordot/Prod:output:0FGN_0/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2
FGN_0/Tensordot/stackе
FGN_0/Tensordot/transpose	Transposeinputs_2FGN_0/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2
FGN_0/Tensordot/transposeи
FGN_0/Tensordot/ReshapeReshapeFGN_0/Tensordot/transpose:y:0FGN_0/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2
FGN_0/Tensordot/ReshapeХ
FGN_0/Tensordot/MatMulMatMul FGN_0/Tensordot/Reshape:output:0&FGN_0/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
FGN_0/Tensordot/MatMul|
FGN_0/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2
FGN_0/Tensordot/Const_2ђ
FGN_0/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_0/Tensordot/concat_1/axis█
FGN_0/Tensordot/concat_1ConcatV2!FGN_0/Tensordot/GatherV2:output:0 FGN_0/Tensordot/Const_2:output:0&FGN_0/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_0/Tensordot/concat_1г
FGN_0/TensordotReshape FGN_0/Tensordot/MatMul:product:0!FGN_0/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
FGN_0/Tensordotъ
FGN_0/BiasAdd/ReadVariableOpReadVariableOp%fgn_0_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
FGN_0/BiasAdd/ReadVariableOpБ
FGN_0/BiasAddBiasAddFGN_0/Tensordot:output:0$FGN_0/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
FGN_0/BiasAddr

FGN_0/ReluReluFGN_0/BiasAdd:output:0*
T0*/
_output_shapes
:         2

FGN_0/Reluе
FGN_1/Tensordot/ReadVariableOpReadVariableOp'fgn_1_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02 
FGN_1/Tensordot/ReadVariableOpv
FGN_1/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2
FGN_1/Tensordot/axesЂ
FGN_1/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2
FGN_1/Tensordot/freev
FGN_1/Tensordot/ShapeShapeFGN_0/Relu:activations:0*
T0*
_output_shapes
:2
FGN_1/Tensordot/Shapeђ
FGN_1/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_1/Tensordot/GatherV2/axis№
FGN_1/Tensordot/GatherV2GatherV2FGN_1/Tensordot/Shape:output:0FGN_1/Tensordot/free:output:0&FGN_1/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_1/Tensordot/GatherV2ё
FGN_1/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2!
FGN_1/Tensordot/GatherV2_1/axisш
FGN_1/Tensordot/GatherV2_1GatherV2FGN_1/Tensordot/Shape:output:0FGN_1/Tensordot/axes:output:0(FGN_1/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_1/Tensordot/GatherV2_1x
FGN_1/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
FGN_1/Tensordot/Constў
FGN_1/Tensordot/ProdProd!FGN_1/Tensordot/GatherV2:output:0FGN_1/Tensordot/Const:output:0*
T0*
_output_shapes
: 2
FGN_1/Tensordot/Prod|
FGN_1/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2
FGN_1/Tensordot/Const_1а
FGN_1/Tensordot/Prod_1Prod#FGN_1/Tensordot/GatherV2_1:output:0 FGN_1/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2
FGN_1/Tensordot/Prod_1|
FGN_1/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_1/Tensordot/concat/axis╬
FGN_1/Tensordot/concatConcatV2FGN_1/Tensordot/free:output:0FGN_1/Tensordot/axes:output:0$FGN_1/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_1/Tensordot/concatц
FGN_1/Tensordot/stackPackFGN_1/Tensordot/Prod:output:0FGN_1/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2
FGN_1/Tensordot/stackИ
FGN_1/Tensordot/transpose	TransposeFGN_0/Relu:activations:0FGN_1/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2
FGN_1/Tensordot/transposeи
FGN_1/Tensordot/ReshapeReshapeFGN_1/Tensordot/transpose:y:0FGN_1/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2
FGN_1/Tensordot/ReshapeХ
FGN_1/Tensordot/MatMulMatMul FGN_1/Tensordot/Reshape:output:0&FGN_1/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
FGN_1/Tensordot/MatMul|
FGN_1/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2
FGN_1/Tensordot/Const_2ђ
FGN_1/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_1/Tensordot/concat_1/axis█
FGN_1/Tensordot/concat_1ConcatV2!FGN_1/Tensordot/GatherV2:output:0 FGN_1/Tensordot/Const_2:output:0&FGN_1/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_1/Tensordot/concat_1г
FGN_1/TensordotReshape FGN_1/Tensordot/MatMul:product:0!FGN_1/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
FGN_1/Tensordotъ
FGN_1/BiasAdd/ReadVariableOpReadVariableOp%fgn_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
FGN_1/BiasAdd/ReadVariableOpБ
FGN_1/BiasAddBiasAddFGN_1/Tensordot:output:0$FGN_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
FGN_1/BiasAddr

FGN_1/ReluReluFGN_1/BiasAdd:output:0*
T0*/
_output_shapes
:         2

FGN_1/Relu«
 FGN_out/Tensordot/ReadVariableOpReadVariableOp)fgn_out_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02"
 FGN_out/Tensordot/ReadVariableOpz
FGN_out/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2
FGN_out/Tensordot/axesЁ
FGN_out/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2
FGN_out/Tensordot/freez
FGN_out/Tensordot/ShapeShapeFGN_1/Relu:activations:0*
T0*
_output_shapes
:2
FGN_out/Tensordot/Shapeё
FGN_out/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2!
FGN_out/Tensordot/GatherV2/axisщ
FGN_out/Tensordot/GatherV2GatherV2 FGN_out/Tensordot/Shape:output:0FGN_out/Tensordot/free:output:0(FGN_out/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_out/Tensordot/GatherV2ѕ
!FGN_out/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2#
!FGN_out/Tensordot/GatherV2_1/axis 
FGN_out/Tensordot/GatherV2_1GatherV2 FGN_out/Tensordot/Shape:output:0FGN_out/Tensordot/axes:output:0*FGN_out/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_out/Tensordot/GatherV2_1|
FGN_out/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
FGN_out/Tensordot/Constа
FGN_out/Tensordot/ProdProd#FGN_out/Tensordot/GatherV2:output:0 FGN_out/Tensordot/Const:output:0*
T0*
_output_shapes
: 2
FGN_out/Tensordot/Prodђ
FGN_out/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2
FGN_out/Tensordot/Const_1е
FGN_out/Tensordot/Prod_1Prod%FGN_out/Tensordot/GatherV2_1:output:0"FGN_out/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2
FGN_out/Tensordot/Prod_1ђ
FGN_out/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_out/Tensordot/concat/axisп
FGN_out/Tensordot/concatConcatV2FGN_out/Tensordot/free:output:0FGN_out/Tensordot/axes:output:0&FGN_out/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_out/Tensordot/concatг
FGN_out/Tensordot/stackPackFGN_out/Tensordot/Prod:output:0!FGN_out/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2
FGN_out/Tensordot/stackЙ
FGN_out/Tensordot/transpose	TransposeFGN_1/Relu:activations:0!FGN_out/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2
FGN_out/Tensordot/transpose┐
FGN_out/Tensordot/ReshapeReshapeFGN_out/Tensordot/transpose:y:0 FGN_out/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2
FGN_out/Tensordot/ReshapeЙ
FGN_out/Tensordot/MatMulMatMul"FGN_out/Tensordot/Reshape:output:0(FGN_out/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
FGN_out/Tensordot/MatMulђ
FGN_out/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2
FGN_out/Tensordot/Const_2ё
FGN_out/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2!
FGN_out/Tensordot/concat_1/axisт
FGN_out/Tensordot/concat_1ConcatV2#FGN_out/Tensordot/GatherV2:output:0"FGN_out/Tensordot/Const_2:output:0(FGN_out/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_out/Tensordot/concat_1┤
FGN_out/TensordotReshape"FGN_out/Tensordot/MatMul:product:0#FGN_out/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
FGN_out/Tensordotц
FGN_out/BiasAdd/ReadVariableOpReadVariableOp'fgn_out_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
FGN_out/BiasAdd/ReadVariableOpФ
FGN_out/BiasAddBiasAddFGN_out/Tensordot:output:0&FGN_out/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
FGN_out/BiasAddm
Reshape/shape/0Const*
_output_shapes
: *
dtype0*
valueB :
         2
Reshape/shape/0d
Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :2
Reshape/shape/3d
Reshape/shape/4Const*
_output_shapes
: *
dtype0*
value	B :2
Reshape/shape/4м
Reshape/shapePackReshape/shape/0:output:0strided_slice:output:0strided_slice:output:0Reshape/shape/3:output:0Reshape/shape/4:output:0*
N*
T0*
_output_shapes
:2
Reshape/shapeЪ
ReshapeReshapeFGN_out/BiasAdd:output:0Reshape/shape:output:0*
T0*E
_output_shapes3
1:/                           2	
ReshapeЃ
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2
strided_slice_1/stackЄ
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"            2
strided_slice_1/stack_1Є
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2
strided_slice_1/stack_2Ћ
strided_slice_1StridedSliceinputs_1strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*3
_output_shapes!
:         *
ellipsis_mask*
new_axis_mask2
strided_slice_1{
mulMulReshape:output:0strided_slice_1:output:0*
T0*3
_output_shapes!
:         2
mulў
einsum/EinsumEinsummul:z:0inputs_0*
N*
T0*+
_output_shapes
:         *
equationabicf,aif->abc2
einsum/EinsumЇ
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOpђ
MatMulBatchMatMulV2inputs_0MatMul/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
MatMulr
addAddV2einsum/Einsum:output:0MatMul:output:0*
T0*+
_output_shapes
:         2
addї
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp|
BiasAddBiasAddadd:z:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2	
BiasAddY
EluEluBiasAdd:output:0*
T0*+
_output_shapes
:         2
Elui
IdentityIdentityElu:activations:0*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*|
_input_shapesk
i:         :         :         :::::::::U Q
+
_output_shapes
:         
"
_user_specified_name
inputs/0:UQ
+
_output_shapes
:         
"
_user_specified_name
inputs/1:YU
/
_output_shapes
:         
"
_user_specified_name
inputs/2:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: 
­
љ
@__inference_MyModel_layer_call_and_return_conditional_losses_784

inputs
inputs_1
inputs_2
conv1layer_761
conv1layer_763
conv1layer_765
conv1layer_767
conv1layer_769
conv1layer_771
conv1layer_773
conv1layer_775
outputlayer_778
outputlayer_780
identityѕб"Conv1Layer/StatefulPartitionedCallб#OutputLayer/StatefulPartitionedCallщ
"Conv1Layer/StatefulPartitionedCallStatefulPartitionedCallinputsinputs_1inputs_2conv1layer_761conv1layer_763conv1layer_765conv1layer_767conv1layer_769conv1layer_771conv1layer_773conv1layer_775*
Tin
2*
Tout
2*+
_output_shapes
:         **
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_Conv1Layer_layer_call_and_return_conditional_losses_6332$
"Conv1Layer/StatefulPartitionedCallА
#OutputLayer/StatefulPartitionedCallStatefulPartitionedCall+Conv1Layer/StatefulPartitionedCall:output:0outputlayer_778outputlayer_780*
Tin
2*
Tout
2*+
_output_shapes
:         *$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*M
fHRF
D__inference_OutputLayer_layer_call_and_return_conditional_losses_7062%
#OutputLayer/StatefulPartitionedCall¤
IdentityIdentity,OutputLayer/StatefulPartitionedCall:output:0#^Conv1Layer/StatefulPartitionedCall$^OutputLayer/StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         ::::::::::2H
"Conv1Layer/StatefulPartitionedCall"Conv1Layer/StatefulPartitionedCall2J
#OutputLayer/StatefulPartitionedCall#OutputLayer/StatefulPartitionedCall:S O
+
_output_shapes
:         
 
_user_specified_nameinputs:SO
+
_output_shapes
:         
 
_user_specified_nameinputs:WS
/
_output_shapes
:         
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
іђ
М
C__inference_Conv1Layer_layer_call_and_return_conditional_losses_633

inputs
inputs_1
inputs_2+
'fgn_0_tensordot_readvariableop_resource)
%fgn_0_biasadd_readvariableop_resource+
'fgn_1_tensordot_readvariableop_resource)
%fgn_1_biasadd_readvariableop_resource-
)fgn_out_tensordot_readvariableop_resource+
'fgn_out_biasadd_readvariableop_resource"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identityѕD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shape}
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2
strided_slice/stackЂ
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2Р
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_sliceе
FGN_0/Tensordot/ReadVariableOpReadVariableOp'fgn_0_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02 
FGN_0/Tensordot/ReadVariableOpv
FGN_0/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2
FGN_0/Tensordot/axesЂ
FGN_0/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2
FGN_0/Tensordot/freef
FGN_0/Tensordot/ShapeShapeinputs_2*
T0*
_output_shapes
:2
FGN_0/Tensordot/Shapeђ
FGN_0/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_0/Tensordot/GatherV2/axis№
FGN_0/Tensordot/GatherV2GatherV2FGN_0/Tensordot/Shape:output:0FGN_0/Tensordot/free:output:0&FGN_0/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_0/Tensordot/GatherV2ё
FGN_0/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2!
FGN_0/Tensordot/GatherV2_1/axisш
FGN_0/Tensordot/GatherV2_1GatherV2FGN_0/Tensordot/Shape:output:0FGN_0/Tensordot/axes:output:0(FGN_0/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_0/Tensordot/GatherV2_1x
FGN_0/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
FGN_0/Tensordot/Constў
FGN_0/Tensordot/ProdProd!FGN_0/Tensordot/GatherV2:output:0FGN_0/Tensordot/Const:output:0*
T0*
_output_shapes
: 2
FGN_0/Tensordot/Prod|
FGN_0/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2
FGN_0/Tensordot/Const_1а
FGN_0/Tensordot/Prod_1Prod#FGN_0/Tensordot/GatherV2_1:output:0 FGN_0/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2
FGN_0/Tensordot/Prod_1|
FGN_0/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_0/Tensordot/concat/axis╬
FGN_0/Tensordot/concatConcatV2FGN_0/Tensordot/free:output:0FGN_0/Tensordot/axes:output:0$FGN_0/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_0/Tensordot/concatц
FGN_0/Tensordot/stackPackFGN_0/Tensordot/Prod:output:0FGN_0/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2
FGN_0/Tensordot/stackе
FGN_0/Tensordot/transpose	Transposeinputs_2FGN_0/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2
FGN_0/Tensordot/transposeи
FGN_0/Tensordot/ReshapeReshapeFGN_0/Tensordot/transpose:y:0FGN_0/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2
FGN_0/Tensordot/ReshapeХ
FGN_0/Tensordot/MatMulMatMul FGN_0/Tensordot/Reshape:output:0&FGN_0/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
FGN_0/Tensordot/MatMul|
FGN_0/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2
FGN_0/Tensordot/Const_2ђ
FGN_0/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_0/Tensordot/concat_1/axis█
FGN_0/Tensordot/concat_1ConcatV2!FGN_0/Tensordot/GatherV2:output:0 FGN_0/Tensordot/Const_2:output:0&FGN_0/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_0/Tensordot/concat_1г
FGN_0/TensordotReshape FGN_0/Tensordot/MatMul:product:0!FGN_0/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
FGN_0/Tensordotъ
FGN_0/BiasAdd/ReadVariableOpReadVariableOp%fgn_0_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
FGN_0/BiasAdd/ReadVariableOpБ
FGN_0/BiasAddBiasAddFGN_0/Tensordot:output:0$FGN_0/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
FGN_0/BiasAddr

FGN_0/ReluReluFGN_0/BiasAdd:output:0*
T0*/
_output_shapes
:         2

FGN_0/Reluе
FGN_1/Tensordot/ReadVariableOpReadVariableOp'fgn_1_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02 
FGN_1/Tensordot/ReadVariableOpv
FGN_1/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2
FGN_1/Tensordot/axesЂ
FGN_1/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2
FGN_1/Tensordot/freev
FGN_1/Tensordot/ShapeShapeFGN_0/Relu:activations:0*
T0*
_output_shapes
:2
FGN_1/Tensordot/Shapeђ
FGN_1/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_1/Tensordot/GatherV2/axis№
FGN_1/Tensordot/GatherV2GatherV2FGN_1/Tensordot/Shape:output:0FGN_1/Tensordot/free:output:0&FGN_1/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_1/Tensordot/GatherV2ё
FGN_1/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2!
FGN_1/Tensordot/GatherV2_1/axisш
FGN_1/Tensordot/GatherV2_1GatherV2FGN_1/Tensordot/Shape:output:0FGN_1/Tensordot/axes:output:0(FGN_1/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_1/Tensordot/GatherV2_1x
FGN_1/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
FGN_1/Tensordot/Constў
FGN_1/Tensordot/ProdProd!FGN_1/Tensordot/GatherV2:output:0FGN_1/Tensordot/Const:output:0*
T0*
_output_shapes
: 2
FGN_1/Tensordot/Prod|
FGN_1/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2
FGN_1/Tensordot/Const_1а
FGN_1/Tensordot/Prod_1Prod#FGN_1/Tensordot/GatherV2_1:output:0 FGN_1/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2
FGN_1/Tensordot/Prod_1|
FGN_1/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_1/Tensordot/concat/axis╬
FGN_1/Tensordot/concatConcatV2FGN_1/Tensordot/free:output:0FGN_1/Tensordot/axes:output:0$FGN_1/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_1/Tensordot/concatц
FGN_1/Tensordot/stackPackFGN_1/Tensordot/Prod:output:0FGN_1/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2
FGN_1/Tensordot/stackИ
FGN_1/Tensordot/transpose	TransposeFGN_0/Relu:activations:0FGN_1/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2
FGN_1/Tensordot/transposeи
FGN_1/Tensordot/ReshapeReshapeFGN_1/Tensordot/transpose:y:0FGN_1/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2
FGN_1/Tensordot/ReshapeХ
FGN_1/Tensordot/MatMulMatMul FGN_1/Tensordot/Reshape:output:0&FGN_1/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
FGN_1/Tensordot/MatMul|
FGN_1/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2
FGN_1/Tensordot/Const_2ђ
FGN_1/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_1/Tensordot/concat_1/axis█
FGN_1/Tensordot/concat_1ConcatV2!FGN_1/Tensordot/GatherV2:output:0 FGN_1/Tensordot/Const_2:output:0&FGN_1/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_1/Tensordot/concat_1г
FGN_1/TensordotReshape FGN_1/Tensordot/MatMul:product:0!FGN_1/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
FGN_1/Tensordotъ
FGN_1/BiasAdd/ReadVariableOpReadVariableOp%fgn_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
FGN_1/BiasAdd/ReadVariableOpБ
FGN_1/BiasAddBiasAddFGN_1/Tensordot:output:0$FGN_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
FGN_1/BiasAddr

FGN_1/ReluReluFGN_1/BiasAdd:output:0*
T0*/
_output_shapes
:         2

FGN_1/Relu«
 FGN_out/Tensordot/ReadVariableOpReadVariableOp)fgn_out_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02"
 FGN_out/Tensordot/ReadVariableOpz
FGN_out/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2
FGN_out/Tensordot/axesЁ
FGN_out/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2
FGN_out/Tensordot/freez
FGN_out/Tensordot/ShapeShapeFGN_1/Relu:activations:0*
T0*
_output_shapes
:2
FGN_out/Tensordot/Shapeё
FGN_out/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2!
FGN_out/Tensordot/GatherV2/axisщ
FGN_out/Tensordot/GatherV2GatherV2 FGN_out/Tensordot/Shape:output:0FGN_out/Tensordot/free:output:0(FGN_out/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_out/Tensordot/GatherV2ѕ
!FGN_out/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2#
!FGN_out/Tensordot/GatherV2_1/axis 
FGN_out/Tensordot/GatherV2_1GatherV2 FGN_out/Tensordot/Shape:output:0FGN_out/Tensordot/axes:output:0*FGN_out/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
FGN_out/Tensordot/GatherV2_1|
FGN_out/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
FGN_out/Tensordot/Constа
FGN_out/Tensordot/ProdProd#FGN_out/Tensordot/GatherV2:output:0 FGN_out/Tensordot/Const:output:0*
T0*
_output_shapes
: 2
FGN_out/Tensordot/Prodђ
FGN_out/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2
FGN_out/Tensordot/Const_1е
FGN_out/Tensordot/Prod_1Prod%FGN_out/Tensordot/GatherV2_1:output:0"FGN_out/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2
FGN_out/Tensordot/Prod_1ђ
FGN_out/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
FGN_out/Tensordot/concat/axisп
FGN_out/Tensordot/concatConcatV2FGN_out/Tensordot/free:output:0FGN_out/Tensordot/axes:output:0&FGN_out/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_out/Tensordot/concatг
FGN_out/Tensordot/stackPackFGN_out/Tensordot/Prod:output:0!FGN_out/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2
FGN_out/Tensordot/stackЙ
FGN_out/Tensordot/transpose	TransposeFGN_1/Relu:activations:0!FGN_out/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2
FGN_out/Tensordot/transpose┐
FGN_out/Tensordot/ReshapeReshapeFGN_out/Tensordot/transpose:y:0 FGN_out/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2
FGN_out/Tensordot/ReshapeЙ
FGN_out/Tensordot/MatMulMatMul"FGN_out/Tensordot/Reshape:output:0(FGN_out/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
FGN_out/Tensordot/MatMulђ
FGN_out/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2
FGN_out/Tensordot/Const_2ё
FGN_out/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2!
FGN_out/Tensordot/concat_1/axisт
FGN_out/Tensordot/concat_1ConcatV2#FGN_out/Tensordot/GatherV2:output:0"FGN_out/Tensordot/Const_2:output:0(FGN_out/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2
FGN_out/Tensordot/concat_1┤
FGN_out/TensordotReshape"FGN_out/Tensordot/MatMul:product:0#FGN_out/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2
FGN_out/Tensordotц
FGN_out/BiasAdd/ReadVariableOpReadVariableOp'fgn_out_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
FGN_out/BiasAdd/ReadVariableOpФ
FGN_out/BiasAddBiasAddFGN_out/Tensordot:output:0&FGN_out/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2
FGN_out/BiasAddm
Reshape/shape/0Const*
_output_shapes
: *
dtype0*
valueB :
         2
Reshape/shape/0d
Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :2
Reshape/shape/3d
Reshape/shape/4Const*
_output_shapes
: *
dtype0*
value	B :2
Reshape/shape/4м
Reshape/shapePackReshape/shape/0:output:0strided_slice:output:0strided_slice:output:0Reshape/shape/3:output:0Reshape/shape/4:output:0*
N*
T0*
_output_shapes
:2
Reshape/shapeЪ
ReshapeReshapeFGN_out/BiasAdd:output:0Reshape/shape:output:0*
T0*E
_output_shapes3
1:/                           2	
ReshapeЃ
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2
strided_slice_1/stackЄ
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"            2
strided_slice_1/stack_1Є
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2
strided_slice_1/stack_2Ћ
strided_slice_1StridedSliceinputs_1strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*3
_output_shapes!
:         *
ellipsis_mask*
new_axis_mask2
strided_slice_1{
mulMulReshape:output:0strided_slice_1:output:0*
T0*3
_output_shapes!
:         2
mulќ
einsum/EinsumEinsummul:z:0inputs*
N*
T0*+
_output_shapes
:         *
equationabicf,aif->abc2
einsum/EinsumЇ
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOp~
MatMulBatchMatMulV2inputsMatMul/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
MatMulr
addAddV2einsum/Einsum:output:0MatMul:output:0*
T0*+
_output_shapes
:         2
addї
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp|
BiasAddBiasAddadd:z:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2	
BiasAddY
EluEluBiasAdd:output:0*
T0*+
_output_shapes
:         2
Elui
IdentityIdentityElu:activations:0*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*|
_input_shapesk
i:         :         :         :::::::::S O
+
_output_shapes
:         
 
_user_specified_nameinputs:SO
+
_output_shapes
:         
 
_user_specified_nameinputs:WS
/
_output_shapes
:         
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: 
ў
»
D__inference_OutputLayer_layer_call_and_return_conditional_losses_706

inputs%
!tensordot_readvariableop_resource#
biasadd_readvariableop_resource
identityѕќ
Tensordot/ReadVariableOpReadVariableOp!tensordot_readvariableop_resource*
_output_shapes

:*
dtype02
Tensordot/ReadVariableOpj
Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2
Tensordot/axesq
Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB"       2
Tensordot/freeX
Tensordot/ShapeShapeinputs*
T0*
_output_shapes
:2
Tensordot/Shapet
Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
Tensordot/GatherV2/axisЛ
Tensordot/GatherV2GatherV2Tensordot/Shape:output:0Tensordot/free:output:0 Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
Tensordot/GatherV2x
Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
Tensordot/GatherV2_1/axisО
Tensordot/GatherV2_1GatherV2Tensordot/Shape:output:0Tensordot/axes:output:0"Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
Tensordot/GatherV2_1l
Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Tensordot/Constђ
Tensordot/ProdProdTensordot/GatherV2:output:0Tensordot/Const:output:0*
T0*
_output_shapes
: 2
Tensordot/Prodp
Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2
Tensordot/Const_1ѕ
Tensordot/Prod_1ProdTensordot/GatherV2_1:output:0Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2
Tensordot/Prod_1p
Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
Tensordot/concat/axis░
Tensordot/concatConcatV2Tensordot/free:output:0Tensordot/axes:output:0Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2
Tensordot/concatї
Tensordot/stackPackTensordot/Prod:output:0Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2
Tensordot/stackљ
Tensordot/transpose	TransposeinputsTensordot/concat:output:0*
T0*+
_output_shapes
:         2
Tensordot/transposeЪ
Tensordot/ReshapeReshapeTensordot/transpose:y:0Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2
Tensordot/Reshapeъ
Tensordot/MatMulMatMulTensordot/Reshape:output:0 Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
Tensordot/MatMulp
Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2
Tensordot/Const_2t
Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
Tensordot/concat_1/axisй
Tensordot/concat_1ConcatV2Tensordot/GatherV2:output:0Tensordot/Const_2:output:0 Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2
Tensordot/concat_1љ
	TensordotReshapeTensordot/MatMul:product:0Tensordot/concat_1:output:0*
T0*+
_output_shapes
:         2
	Tensordotї
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOpЄ
BiasAddBiasAddTensordot:output:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2	
BiasAddY
EluEluBiasAdd:output:0*
T0*+
_output_shapes
:         2
Elui
IdentityIdentityElu:activations:0*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*2
_input_shapes!
:         :::S O
+
_output_shapes
:         
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
У
ё
%__inference_MyModel_layer_call_fn_862
x_in
a_in
e_in
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
identityѕбStatefulPartitionedCallк
StatefulPartitionedCallStatefulPartitionedCallx_ina_ine_inunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8*
Tin
2*
Tout
2*+
_output_shapes
:         *,
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8*I
fDRB
@__inference_MyModel_layer_call_and_return_conditional_losses_8392
StatefulPartitionedCallњ
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         ::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:Q M
+
_output_shapes
:         

_user_specified_nameX_in:QM
+
_output_shapes
:         

_user_specified_nameA_in:UQ
/
_output_shapes
:         

_user_specified_nameE_in:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
З2
Э
 __inference__traced_restore_1493
file_prefix+
'assignvariableop_conv1layer_root_kernel&
"assignvariableop_1_conv1layer_bias)
%assignvariableop_2_outputlayer_kernel'
#assignvariableop_3_outputlayer_bias.
*assignvariableop_4_conv1layer_fgn_0_kernel,
(assignvariableop_5_conv1layer_fgn_0_bias.
*assignvariableop_6_conv1layer_fgn_1_kernel,
(assignvariableop_7_conv1layer_fgn_1_bias0
,assignvariableop_8_conv1layer_fgn_out_kernel.
*assignvariableop_9_conv1layer_fgn_out_bias
identity_11ѕбAssignVariableOpбAssignVariableOp_1бAssignVariableOp_2бAssignVariableOp_3бAssignVariableOp_4бAssignVariableOp_5бAssignVariableOp_6бAssignVariableOp_7бAssignVariableOp_8бAssignVariableOp_9б	RestoreV2бRestoreV2_1ќ
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:
*
dtype0*б
valueўBЋ
B;layer_with_weights-0/root_kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUE2
RestoreV2/tensor_namesб
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:
*
dtype0*'
valueB
B B B B B B B B B B 2
RestoreV2/shape_and_slicesП
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*<
_output_shapes*
(::::::::::*
dtypes
2
2
	RestoreV2X
IdentityIdentityRestoreV2:tensors:0*
T0*
_output_shapes
:2

IdentityЌ
AssignVariableOpAssignVariableOp'assignvariableop_conv1layer_root_kernelIdentity:output:0*
_output_shapes
 *
dtype02
AssignVariableOp\

Identity_1IdentityRestoreV2:tensors:1*
T0*
_output_shapes
:2

Identity_1ў
AssignVariableOp_1AssignVariableOp"assignvariableop_1_conv1layer_biasIdentity_1:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_1\

Identity_2IdentityRestoreV2:tensors:2*
T0*
_output_shapes
:2

Identity_2Џ
AssignVariableOp_2AssignVariableOp%assignvariableop_2_outputlayer_kernelIdentity_2:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_2\

Identity_3IdentityRestoreV2:tensors:3*
T0*
_output_shapes
:2

Identity_3Ў
AssignVariableOp_3AssignVariableOp#assignvariableop_3_outputlayer_biasIdentity_3:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_3\

Identity_4IdentityRestoreV2:tensors:4*
T0*
_output_shapes
:2

Identity_4а
AssignVariableOp_4AssignVariableOp*assignvariableop_4_conv1layer_fgn_0_kernelIdentity_4:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_4\

Identity_5IdentityRestoreV2:tensors:5*
T0*
_output_shapes
:2

Identity_5ъ
AssignVariableOp_5AssignVariableOp(assignvariableop_5_conv1layer_fgn_0_biasIdentity_5:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_5\

Identity_6IdentityRestoreV2:tensors:6*
T0*
_output_shapes
:2

Identity_6а
AssignVariableOp_6AssignVariableOp*assignvariableop_6_conv1layer_fgn_1_kernelIdentity_6:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_6\

Identity_7IdentityRestoreV2:tensors:7*
T0*
_output_shapes
:2

Identity_7ъ
AssignVariableOp_7AssignVariableOp(assignvariableop_7_conv1layer_fgn_1_biasIdentity_7:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_7\

Identity_8IdentityRestoreV2:tensors:8*
T0*
_output_shapes
:2

Identity_8б
AssignVariableOp_8AssignVariableOp,assignvariableop_8_conv1layer_fgn_out_kernelIdentity_8:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_8\

Identity_9IdentityRestoreV2:tensors:9*
T0*
_output_shapes
:2

Identity_9а
AssignVariableOp_9AssignVariableOp*assignvariableop_9_conv1layer_fgn_out_biasIdentity_9:output:0*
_output_shapes
 *
dtype02
AssignVariableOp_9е
RestoreV2_1/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*1
value(B&B_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2_1/tensor_namesћ
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*
valueB
B 2
RestoreV2_1/shape_and_slices─
RestoreV2_1	RestoreV2file_prefix!RestoreV2_1/tensor_names:output:0%RestoreV2_1/shape_and_slices:output:0
^RestoreV2"/device:CPU:0*
_output_shapes
:*
dtypes
22
RestoreV2_19
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp║
Identity_10Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_10К
Identity_11IdentityIdentity_10:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9
^RestoreV2^RestoreV2_1*
T0*
_output_shapes
: 2
Identity_11"#
identity_11Identity_11:output:0*=
_input_shapes,
*: ::::::::::2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12(
AssignVariableOp_2AssignVariableOp_22(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_92
	RestoreV2	RestoreV22
RestoreV2_1RestoreV2_1:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: 
Ї
Љ
&__inference_MyModel_layer_call_fn_1219
inputs_0
inputs_1
inputs_2
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
identityѕбStatefulPartitionedCallм
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1inputs_2unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8*
Tin
2*
Tout
2*+
_output_shapes
:         *,
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8*I
fDRB
@__inference_MyModel_layer_call_and_return_conditional_losses_8392
StatefulPartitionedCallњ
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         ::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:U Q
+
_output_shapes
:         
"
_user_specified_name
inputs/0:UQ
+
_output_shapes
:         
"
_user_specified_name
inputs/1:YU
/
_output_shapes
:         
"
_user_specified_name
inputs/2:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
о
є
@__inference_MyModel_layer_call_and_return_conditional_losses_723
x_in
a_in
e_in
conv1layer_658
conv1layer_660
conv1layer_662
conv1layer_664
conv1layer_666
conv1layer_668
conv1layer_670
conv1layer_672
outputlayer_717
outputlayer_719
identityѕб"Conv1Layer/StatefulPartitionedCallб#OutputLayer/StatefulPartitionedCall№
"Conv1Layer/StatefulPartitionedCallStatefulPartitionedCallx_ina_ine_inconv1layer_658conv1layer_660conv1layer_662conv1layer_664conv1layer_666conv1layer_668conv1layer_670conv1layer_672*
Tin
2*
Tout
2*+
_output_shapes
:         **
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_Conv1Layer_layer_call_and_return_conditional_losses_6332$
"Conv1Layer/StatefulPartitionedCallА
#OutputLayer/StatefulPartitionedCallStatefulPartitionedCall+Conv1Layer/StatefulPartitionedCall:output:0outputlayer_717outputlayer_719*
Tin
2*
Tout
2*+
_output_shapes
:         *$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*M
fHRF
D__inference_OutputLayer_layer_call_and_return_conditional_losses_7062%
#OutputLayer/StatefulPartitionedCall¤
IdentityIdentity,OutputLayer/StatefulPartitionedCall:output:0#^Conv1Layer/StatefulPartitionedCall$^OutputLayer/StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         ::::::::::2H
"Conv1Layer/StatefulPartitionedCall"Conv1Layer/StatefulPartitionedCall2J
#OutputLayer/StatefulPartitionedCall#OutputLayer/StatefulPartitionedCall:Q M
+
_output_shapes
:         

_user_specified_nameX_in:QM
+
_output_shapes
:         

_user_specified_nameA_in:UQ
/
_output_shapes
:         

_user_specified_nameE_in:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
о
є
@__inference_MyModel_layer_call_and_return_conditional_losses_751
x_in
a_in
e_in
conv1layer_728
conv1layer_730
conv1layer_732
conv1layer_734
conv1layer_736
conv1layer_738
conv1layer_740
conv1layer_742
outputlayer_745
outputlayer_747
identityѕб"Conv1Layer/StatefulPartitionedCallб#OutputLayer/StatefulPartitionedCall№
"Conv1Layer/StatefulPartitionedCallStatefulPartitionedCallx_ina_ine_inconv1layer_728conv1layer_730conv1layer_732conv1layer_734conv1layer_736conv1layer_738conv1layer_740conv1layer_742*
Tin
2*
Tout
2*+
_output_shapes
:         **
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_Conv1Layer_layer_call_and_return_conditional_losses_6332$
"Conv1Layer/StatefulPartitionedCallА
#OutputLayer/StatefulPartitionedCallStatefulPartitionedCall+Conv1Layer/StatefulPartitionedCall:output:0outputlayer_745outputlayer_747*
Tin
2*
Tout
2*+
_output_shapes
:         *$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*M
fHRF
D__inference_OutputLayer_layer_call_and_return_conditional_losses_7062%
#OutputLayer/StatefulPartitionedCall¤
IdentityIdentity,OutputLayer/StatefulPartitionedCall:output:0#^Conv1Layer/StatefulPartitionedCall$^OutputLayer/StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         ::::::::::2H
"Conv1Layer/StatefulPartitionedCall"Conv1Layer/StatefulPartitionedCall2J
#OutputLayer/StatefulPartitionedCall#OutputLayer/StatefulPartitionedCall:Q M
+
_output_shapes
:         

_user_specified_nameX_in:QM
+
_output_shapes
:         

_user_specified_nameA_in:UQ
/
_output_shapes
:         

_user_specified_nameE_in:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
┬
ђ
!__inference_signature_wrapper_891
a_in
e_in
x_in
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
identityѕбStatefulPartitionedCallц
StatefulPartitionedCallStatefulPartitionedCallx_ina_ine_inunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8*
Tin
2*
Tout
2*+
_output_shapes
:         *,
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8*'
f"R 
__inference__wrapped_model_5172
StatefulPartitionedCallњ
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         ::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:Q M
+
_output_shapes
:         

_user_specified_nameA_in:UQ
/
_output_shapes
:         

_user_specified_nameE_in:QM
+
_output_shapes
:         

_user_specified_nameX_in:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
­
љ
@__inference_MyModel_layer_call_and_return_conditional_losses_839

inputs
inputs_1
inputs_2
conv1layer_816
conv1layer_818
conv1layer_820
conv1layer_822
conv1layer_824
conv1layer_826
conv1layer_828
conv1layer_830
outputlayer_833
outputlayer_835
identityѕб"Conv1Layer/StatefulPartitionedCallб#OutputLayer/StatefulPartitionedCallщ
"Conv1Layer/StatefulPartitionedCallStatefulPartitionedCallinputsinputs_1inputs_2conv1layer_816conv1layer_818conv1layer_820conv1layer_822conv1layer_824conv1layer_826conv1layer_828conv1layer_830*
Tin
2*
Tout
2*+
_output_shapes
:         **
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_Conv1Layer_layer_call_and_return_conditional_losses_6332$
"Conv1Layer/StatefulPartitionedCallА
#OutputLayer/StatefulPartitionedCallStatefulPartitionedCall+Conv1Layer/StatefulPartitionedCall:output:0outputlayer_833outputlayer_835*
Tin
2*
Tout
2*+
_output_shapes
:         *$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*M
fHRF
D__inference_OutputLayer_layer_call_and_return_conditional_losses_7062%
#OutputLayer/StatefulPartitionedCall¤
IdentityIdentity,OutputLayer/StatefulPartitionedCall:output:0#^Conv1Layer/StatefulPartitionedCall$^OutputLayer/StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         ::::::::::2H
"Conv1Layer/StatefulPartitionedCall"Conv1Layer/StatefulPartitionedCall2J
#OutputLayer/StatefulPartitionedCall#OutputLayer/StatefulPartitionedCall:S O
+
_output_shapes
:         
 
_user_specified_nameinputs:SO
+
_output_shapes
:         
 
_user_specified_nameinputs:WS
/
_output_shapes
:         
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
Ў
░
E__inference_OutputLayer_layer_call_and_return_conditional_losses_1383

inputs%
!tensordot_readvariableop_resource#
biasadd_readvariableop_resource
identityѕќ
Tensordot/ReadVariableOpReadVariableOp!tensordot_readvariableop_resource*
_output_shapes

:*
dtype02
Tensordot/ReadVariableOpj
Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2
Tensordot/axesq
Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB"       2
Tensordot/freeX
Tensordot/ShapeShapeinputs*
T0*
_output_shapes
:2
Tensordot/Shapet
Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
Tensordot/GatherV2/axisЛ
Tensordot/GatherV2GatherV2Tensordot/Shape:output:0Tensordot/free:output:0 Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
Tensordot/GatherV2x
Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
Tensordot/GatherV2_1/axisО
Tensordot/GatherV2_1GatherV2Tensordot/Shape:output:0Tensordot/axes:output:0"Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2
Tensordot/GatherV2_1l
Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2
Tensordot/Constђ
Tensordot/ProdProdTensordot/GatherV2:output:0Tensordot/Const:output:0*
T0*
_output_shapes
: 2
Tensordot/Prodp
Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2
Tensordot/Const_1ѕ
Tensordot/Prod_1ProdTensordot/GatherV2_1:output:0Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2
Tensordot/Prod_1p
Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
Tensordot/concat/axis░
Tensordot/concatConcatV2Tensordot/free:output:0Tensordot/axes:output:0Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2
Tensordot/concatї
Tensordot/stackPackTensordot/Prod:output:0Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2
Tensordot/stackљ
Tensordot/transpose	TransposeinputsTensordot/concat:output:0*
T0*+
_output_shapes
:         2
Tensordot/transposeЪ
Tensordot/ReshapeReshapeTensordot/transpose:y:0Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2
Tensordot/Reshapeъ
Tensordot/MatMulMatMulTensordot/Reshape:output:0 Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2
Tensordot/MatMulp
Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2
Tensordot/Const_2t
Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2
Tensordot/concat_1/axisй
Tensordot/concat_1ConcatV2Tensordot/GatherV2:output:0Tensordot/Const_2:output:0 Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2
Tensordot/concat_1љ
	TensordotReshapeTensordot/MatMul:product:0Tensordot/concat_1:output:0*
T0*+
_output_shapes
:         2
	Tensordotї
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOpЄ
BiasAddBiasAddTensordot:output:0BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2	
BiasAddY
EluEluBiasAdd:output:0*
T0*+
_output_shapes
:         2
Elui
IdentityIdentityElu:activations:0*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*2
_input_shapes!
:         :::S O
+
_output_shapes
:         
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
Є

*__inference_OutputLayer_layer_call_fn_1392

inputs
unknown
	unknown_0
identityѕбStatefulPartitionedCallо
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*+
_output_shapes
:         *$
_read_only_resource_inputs
**
config_proto

CPU

GPU 2J 8*M
fHRF
D__inference_OutputLayer_layer_call_and_return_conditional_losses_7062
StatefulPartitionedCallњ
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*2
_input_shapes!
:         ::22
StatefulPartitionedCallStatefulPartitionedCall:S O
+
_output_shapes
:         
 
_user_specified_nameinputs:

_output_shapes
: :

_output_shapes
: 
У
ё
%__inference_MyModel_layer_call_fn_807
x_in
a_in
e_in
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
identityѕбStatefulPartitionedCallк
StatefulPartitionedCallStatefulPartitionedCallx_ina_ine_inunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8*
Tin
2*
Tout
2*+
_output_shapes
:         *,
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8*I
fDRB
@__inference_MyModel_layer_call_and_return_conditional_losses_7842
StatefulPartitionedCallњ
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         ::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:Q M
+
_output_shapes
:         

_user_specified_nameX_in:QM
+
_output_shapes
:         

_user_specified_nameA_in:UQ
/
_output_shapes
:         

_user_specified_nameE_in:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
сп
░
__inference__wrapped_model_517
x_in
a_in
e_in>
:mymodel_conv1layer_fgn_0_tensordot_readvariableop_resource<
8mymodel_conv1layer_fgn_0_biasadd_readvariableop_resource>
:mymodel_conv1layer_fgn_1_tensordot_readvariableop_resource<
8mymodel_conv1layer_fgn_1_biasadd_readvariableop_resource@
<mymodel_conv1layer_fgn_out_tensordot_readvariableop_resource>
:mymodel_conv1layer_fgn_out_biasadd_readvariableop_resource5
1mymodel_conv1layer_matmul_readvariableop_resource6
2mymodel_conv1layer_biasadd_readvariableop_resource9
5mymodel_outputlayer_tensordot_readvariableop_resource7
3mymodel_outputlayer_biasadd_readvariableop_resource
identityѕh
MyModel/Conv1Layer/ShapeShapex_in*
T0*
_output_shapes
:2
MyModel/Conv1Layer/ShapeБ
&MyModel/Conv1Layer/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:
■        2(
&MyModel/Conv1Layer/strided_slice/stackД
(MyModel/Conv1Layer/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:
         2*
(MyModel/Conv1Layer/strided_slice/stack_1ъ
(MyModel/Conv1Layer/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(MyModel/Conv1Layer/strided_slice/stack_2н
 MyModel/Conv1Layer/strided_sliceStridedSlice!MyModel/Conv1Layer/Shape:output:0/MyModel/Conv1Layer/strided_slice/stack:output:01MyModel/Conv1Layer/strided_slice/stack_1:output:01MyModel/Conv1Layer/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2"
 MyModel/Conv1Layer/strided_sliceр
1MyModel/Conv1Layer/FGN_0/Tensordot/ReadVariableOpReadVariableOp:mymodel_conv1layer_fgn_0_tensordot_readvariableop_resource*
_output_shapes

:*
dtype023
1MyModel/Conv1Layer/FGN_0/Tensordot/ReadVariableOpю
'MyModel/Conv1Layer/FGN_0/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2)
'MyModel/Conv1Layer/FGN_0/Tensordot/axesД
'MyModel/Conv1Layer/FGN_0/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2)
'MyModel/Conv1Layer/FGN_0/Tensordot/freeѕ
(MyModel/Conv1Layer/FGN_0/Tensordot/ShapeShapee_in*
T0*
_output_shapes
:2*
(MyModel/Conv1Layer/FGN_0/Tensordot/Shapeд
0MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 22
0MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2/axis╬
+MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2GatherV21MyModel/Conv1Layer/FGN_0/Tensordot/Shape:output:00MyModel/Conv1Layer/FGN_0/Tensordot/free:output:09MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2-
+MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2ф
2MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 24
2MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2_1/axisн
-MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2_1GatherV21MyModel/Conv1Layer/FGN_0/Tensordot/Shape:output:00MyModel/Conv1Layer/FGN_0/Tensordot/axes:output:0;MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2/
-MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2_1ъ
(MyModel/Conv1Layer/FGN_0/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2*
(MyModel/Conv1Layer/FGN_0/Tensordot/ConstС
'MyModel/Conv1Layer/FGN_0/Tensordot/ProdProd4MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2:output:01MyModel/Conv1Layer/FGN_0/Tensordot/Const:output:0*
T0*
_output_shapes
: 2)
'MyModel/Conv1Layer/FGN_0/Tensordot/Prodб
*MyModel/Conv1Layer/FGN_0/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2,
*MyModel/Conv1Layer/FGN_0/Tensordot/Const_1В
)MyModel/Conv1Layer/FGN_0/Tensordot/Prod_1Prod6MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2_1:output:03MyModel/Conv1Layer/FGN_0/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2+
)MyModel/Conv1Layer/FGN_0/Tensordot/Prod_1б
.MyModel/Conv1Layer/FGN_0/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 20
.MyModel/Conv1Layer/FGN_0/Tensordot/concat/axisГ
)MyModel/Conv1Layer/FGN_0/Tensordot/concatConcatV20MyModel/Conv1Layer/FGN_0/Tensordot/free:output:00MyModel/Conv1Layer/FGN_0/Tensordot/axes:output:07MyModel/Conv1Layer/FGN_0/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2+
)MyModel/Conv1Layer/FGN_0/Tensordot/concat­
(MyModel/Conv1Layer/FGN_0/Tensordot/stackPack0MyModel/Conv1Layer/FGN_0/Tensordot/Prod:output:02MyModel/Conv1Layer/FGN_0/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2*
(MyModel/Conv1Layer/FGN_0/Tensordot/stackП
,MyModel/Conv1Layer/FGN_0/Tensordot/transpose	Transposee_in2MyModel/Conv1Layer/FGN_0/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2.
,MyModel/Conv1Layer/FGN_0/Tensordot/transposeЃ
*MyModel/Conv1Layer/FGN_0/Tensordot/ReshapeReshape0MyModel/Conv1Layer/FGN_0/Tensordot/transpose:y:01MyModel/Conv1Layer/FGN_0/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2,
*MyModel/Conv1Layer/FGN_0/Tensordot/Reshapeѓ
)MyModel/Conv1Layer/FGN_0/Tensordot/MatMulMatMul3MyModel/Conv1Layer/FGN_0/Tensordot/Reshape:output:09MyModel/Conv1Layer/FGN_0/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2+
)MyModel/Conv1Layer/FGN_0/Tensordot/MatMulб
*MyModel/Conv1Layer/FGN_0/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*MyModel/Conv1Layer/FGN_0/Tensordot/Const_2д
0MyModel/Conv1Layer/FGN_0/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 22
0MyModel/Conv1Layer/FGN_0/Tensordot/concat_1/axis║
+MyModel/Conv1Layer/FGN_0/Tensordot/concat_1ConcatV24MyModel/Conv1Layer/FGN_0/Tensordot/GatherV2:output:03MyModel/Conv1Layer/FGN_0/Tensordot/Const_2:output:09MyModel/Conv1Layer/FGN_0/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2-
+MyModel/Conv1Layer/FGN_0/Tensordot/concat_1Э
"MyModel/Conv1Layer/FGN_0/TensordotReshape3MyModel/Conv1Layer/FGN_0/Tensordot/MatMul:product:04MyModel/Conv1Layer/FGN_0/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2$
"MyModel/Conv1Layer/FGN_0/TensordotО
/MyModel/Conv1Layer/FGN_0/BiasAdd/ReadVariableOpReadVariableOp8mymodel_conv1layer_fgn_0_biasadd_readvariableop_resource*
_output_shapes
:*
dtype021
/MyModel/Conv1Layer/FGN_0/BiasAdd/ReadVariableOp№
 MyModel/Conv1Layer/FGN_0/BiasAddBiasAdd+MyModel/Conv1Layer/FGN_0/Tensordot:output:07MyModel/Conv1Layer/FGN_0/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2"
 MyModel/Conv1Layer/FGN_0/BiasAddФ
MyModel/Conv1Layer/FGN_0/ReluRelu)MyModel/Conv1Layer/FGN_0/BiasAdd:output:0*
T0*/
_output_shapes
:         2
MyModel/Conv1Layer/FGN_0/Reluр
1MyModel/Conv1Layer/FGN_1/Tensordot/ReadVariableOpReadVariableOp:mymodel_conv1layer_fgn_1_tensordot_readvariableop_resource*
_output_shapes

:*
dtype023
1MyModel/Conv1Layer/FGN_1/Tensordot/ReadVariableOpю
'MyModel/Conv1Layer/FGN_1/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2)
'MyModel/Conv1Layer/FGN_1/Tensordot/axesД
'MyModel/Conv1Layer/FGN_1/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2)
'MyModel/Conv1Layer/FGN_1/Tensordot/free»
(MyModel/Conv1Layer/FGN_1/Tensordot/ShapeShape+MyModel/Conv1Layer/FGN_0/Relu:activations:0*
T0*
_output_shapes
:2*
(MyModel/Conv1Layer/FGN_1/Tensordot/Shapeд
0MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 22
0MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2/axis╬
+MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2GatherV21MyModel/Conv1Layer/FGN_1/Tensordot/Shape:output:00MyModel/Conv1Layer/FGN_1/Tensordot/free:output:09MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2-
+MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2ф
2MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 24
2MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2_1/axisн
-MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2_1GatherV21MyModel/Conv1Layer/FGN_1/Tensordot/Shape:output:00MyModel/Conv1Layer/FGN_1/Tensordot/axes:output:0;MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2/
-MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2_1ъ
(MyModel/Conv1Layer/FGN_1/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2*
(MyModel/Conv1Layer/FGN_1/Tensordot/ConstС
'MyModel/Conv1Layer/FGN_1/Tensordot/ProdProd4MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2:output:01MyModel/Conv1Layer/FGN_1/Tensordot/Const:output:0*
T0*
_output_shapes
: 2)
'MyModel/Conv1Layer/FGN_1/Tensordot/Prodб
*MyModel/Conv1Layer/FGN_1/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2,
*MyModel/Conv1Layer/FGN_1/Tensordot/Const_1В
)MyModel/Conv1Layer/FGN_1/Tensordot/Prod_1Prod6MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2_1:output:03MyModel/Conv1Layer/FGN_1/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2+
)MyModel/Conv1Layer/FGN_1/Tensordot/Prod_1б
.MyModel/Conv1Layer/FGN_1/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 20
.MyModel/Conv1Layer/FGN_1/Tensordot/concat/axisГ
)MyModel/Conv1Layer/FGN_1/Tensordot/concatConcatV20MyModel/Conv1Layer/FGN_1/Tensordot/free:output:00MyModel/Conv1Layer/FGN_1/Tensordot/axes:output:07MyModel/Conv1Layer/FGN_1/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2+
)MyModel/Conv1Layer/FGN_1/Tensordot/concat­
(MyModel/Conv1Layer/FGN_1/Tensordot/stackPack0MyModel/Conv1Layer/FGN_1/Tensordot/Prod:output:02MyModel/Conv1Layer/FGN_1/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2*
(MyModel/Conv1Layer/FGN_1/Tensordot/stackё
,MyModel/Conv1Layer/FGN_1/Tensordot/transpose	Transpose+MyModel/Conv1Layer/FGN_0/Relu:activations:02MyModel/Conv1Layer/FGN_1/Tensordot/concat:output:0*
T0*/
_output_shapes
:         2.
,MyModel/Conv1Layer/FGN_1/Tensordot/transposeЃ
*MyModel/Conv1Layer/FGN_1/Tensordot/ReshapeReshape0MyModel/Conv1Layer/FGN_1/Tensordot/transpose:y:01MyModel/Conv1Layer/FGN_1/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2,
*MyModel/Conv1Layer/FGN_1/Tensordot/Reshapeѓ
)MyModel/Conv1Layer/FGN_1/Tensordot/MatMulMatMul3MyModel/Conv1Layer/FGN_1/Tensordot/Reshape:output:09MyModel/Conv1Layer/FGN_1/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2+
)MyModel/Conv1Layer/FGN_1/Tensordot/MatMulб
*MyModel/Conv1Layer/FGN_1/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*MyModel/Conv1Layer/FGN_1/Tensordot/Const_2д
0MyModel/Conv1Layer/FGN_1/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 22
0MyModel/Conv1Layer/FGN_1/Tensordot/concat_1/axis║
+MyModel/Conv1Layer/FGN_1/Tensordot/concat_1ConcatV24MyModel/Conv1Layer/FGN_1/Tensordot/GatherV2:output:03MyModel/Conv1Layer/FGN_1/Tensordot/Const_2:output:09MyModel/Conv1Layer/FGN_1/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2-
+MyModel/Conv1Layer/FGN_1/Tensordot/concat_1Э
"MyModel/Conv1Layer/FGN_1/TensordotReshape3MyModel/Conv1Layer/FGN_1/Tensordot/MatMul:product:04MyModel/Conv1Layer/FGN_1/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2$
"MyModel/Conv1Layer/FGN_1/TensordotО
/MyModel/Conv1Layer/FGN_1/BiasAdd/ReadVariableOpReadVariableOp8mymodel_conv1layer_fgn_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype021
/MyModel/Conv1Layer/FGN_1/BiasAdd/ReadVariableOp№
 MyModel/Conv1Layer/FGN_1/BiasAddBiasAdd+MyModel/Conv1Layer/FGN_1/Tensordot:output:07MyModel/Conv1Layer/FGN_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2"
 MyModel/Conv1Layer/FGN_1/BiasAddФ
MyModel/Conv1Layer/FGN_1/ReluRelu)MyModel/Conv1Layer/FGN_1/BiasAdd:output:0*
T0*/
_output_shapes
:         2
MyModel/Conv1Layer/FGN_1/Reluу
3MyModel/Conv1Layer/FGN_out/Tensordot/ReadVariableOpReadVariableOp<mymodel_conv1layer_fgn_out_tensordot_readvariableop_resource*
_output_shapes

:*
dtype025
3MyModel/Conv1Layer/FGN_out/Tensordot/ReadVariableOpа
)MyModel/Conv1Layer/FGN_out/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2+
)MyModel/Conv1Layer/FGN_out/Tensordot/axesФ
)MyModel/Conv1Layer/FGN_out/Tensordot/freeConst*
_output_shapes
:*
dtype0*!
valueB"          2+
)MyModel/Conv1Layer/FGN_out/Tensordot/free│
*MyModel/Conv1Layer/FGN_out/Tensordot/ShapeShape+MyModel/Conv1Layer/FGN_1/Relu:activations:0*
T0*
_output_shapes
:2,
*MyModel/Conv1Layer/FGN_out/Tensordot/Shapeф
2MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 24
2MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2/axisп
-MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2GatherV23MyModel/Conv1Layer/FGN_out/Tensordot/Shape:output:02MyModel/Conv1Layer/FGN_out/Tensordot/free:output:0;MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2/
-MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2«
4MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 26
4MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2_1/axisя
/MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2_1GatherV23MyModel/Conv1Layer/FGN_out/Tensordot/Shape:output:02MyModel/Conv1Layer/FGN_out/Tensordot/axes:output:0=MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:21
/MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2_1б
*MyModel/Conv1Layer/FGN_out/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2,
*MyModel/Conv1Layer/FGN_out/Tensordot/ConstВ
)MyModel/Conv1Layer/FGN_out/Tensordot/ProdProd6MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2:output:03MyModel/Conv1Layer/FGN_out/Tensordot/Const:output:0*
T0*
_output_shapes
: 2+
)MyModel/Conv1Layer/FGN_out/Tensordot/Prodд
,MyModel/Conv1Layer/FGN_out/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2.
,MyModel/Conv1Layer/FGN_out/Tensordot/Const_1З
+MyModel/Conv1Layer/FGN_out/Tensordot/Prod_1Prod8MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2_1:output:05MyModel/Conv1Layer/FGN_out/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2-
+MyModel/Conv1Layer/FGN_out/Tensordot/Prod_1д
0MyModel/Conv1Layer/FGN_out/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 22
0MyModel/Conv1Layer/FGN_out/Tensordot/concat/axisи
+MyModel/Conv1Layer/FGN_out/Tensordot/concatConcatV22MyModel/Conv1Layer/FGN_out/Tensordot/free:output:02MyModel/Conv1Layer/FGN_out/Tensordot/axes:output:09MyModel/Conv1Layer/FGN_out/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2-
+MyModel/Conv1Layer/FGN_out/Tensordot/concatЭ
*MyModel/Conv1Layer/FGN_out/Tensordot/stackPack2MyModel/Conv1Layer/FGN_out/Tensordot/Prod:output:04MyModel/Conv1Layer/FGN_out/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2,
*MyModel/Conv1Layer/FGN_out/Tensordot/stackі
.MyModel/Conv1Layer/FGN_out/Tensordot/transpose	Transpose+MyModel/Conv1Layer/FGN_1/Relu:activations:04MyModel/Conv1Layer/FGN_out/Tensordot/concat:output:0*
T0*/
_output_shapes
:         20
.MyModel/Conv1Layer/FGN_out/Tensordot/transposeІ
,MyModel/Conv1Layer/FGN_out/Tensordot/ReshapeReshape2MyModel/Conv1Layer/FGN_out/Tensordot/transpose:y:03MyModel/Conv1Layer/FGN_out/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2.
,MyModel/Conv1Layer/FGN_out/Tensordot/Reshapeі
+MyModel/Conv1Layer/FGN_out/Tensordot/MatMulMatMul5MyModel/Conv1Layer/FGN_out/Tensordot/Reshape:output:0;MyModel/Conv1Layer/FGN_out/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2-
+MyModel/Conv1Layer/FGN_out/Tensordot/MatMulд
,MyModel/Conv1Layer/FGN_out/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2.
,MyModel/Conv1Layer/FGN_out/Tensordot/Const_2ф
2MyModel/Conv1Layer/FGN_out/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 24
2MyModel/Conv1Layer/FGN_out/Tensordot/concat_1/axis─
-MyModel/Conv1Layer/FGN_out/Tensordot/concat_1ConcatV26MyModel/Conv1Layer/FGN_out/Tensordot/GatherV2:output:05MyModel/Conv1Layer/FGN_out/Tensordot/Const_2:output:0;MyModel/Conv1Layer/FGN_out/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2/
-MyModel/Conv1Layer/FGN_out/Tensordot/concat_1ђ
$MyModel/Conv1Layer/FGN_out/TensordotReshape5MyModel/Conv1Layer/FGN_out/Tensordot/MatMul:product:06MyModel/Conv1Layer/FGN_out/Tensordot/concat_1:output:0*
T0*/
_output_shapes
:         2&
$MyModel/Conv1Layer/FGN_out/TensordotП
1MyModel/Conv1Layer/FGN_out/BiasAdd/ReadVariableOpReadVariableOp:mymodel_conv1layer_fgn_out_biasadd_readvariableop_resource*
_output_shapes
:*
dtype023
1MyModel/Conv1Layer/FGN_out/BiasAdd/ReadVariableOpэ
"MyModel/Conv1Layer/FGN_out/BiasAddBiasAdd-MyModel/Conv1Layer/FGN_out/Tensordot:output:09MyModel/Conv1Layer/FGN_out/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:         2$
"MyModel/Conv1Layer/FGN_out/BiasAddЊ
"MyModel/Conv1Layer/Reshape/shape/0Const*
_output_shapes
: *
dtype0*
valueB :
         2$
"MyModel/Conv1Layer/Reshape/shape/0і
"MyModel/Conv1Layer/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :2$
"MyModel/Conv1Layer/Reshape/shape/3і
"MyModel/Conv1Layer/Reshape/shape/4Const*
_output_shapes
: *
dtype0*
value	B :2$
"MyModel/Conv1Layer/Reshape/shape/4О
 MyModel/Conv1Layer/Reshape/shapePack+MyModel/Conv1Layer/Reshape/shape/0:output:0)MyModel/Conv1Layer/strided_slice:output:0)MyModel/Conv1Layer/strided_slice:output:0+MyModel/Conv1Layer/Reshape/shape/3:output:0+MyModel/Conv1Layer/Reshape/shape/4:output:0*
N*
T0*
_output_shapes
:2"
 MyModel/Conv1Layer/Reshape/shapeв
MyModel/Conv1Layer/ReshapeReshape+MyModel/Conv1Layer/FGN_out/BiasAdd:output:0)MyModel/Conv1Layer/Reshape/shape:output:0*
T0*E
_output_shapes3
1:/                           2
MyModel/Conv1Layer/ReshapeЕ
(MyModel/Conv1Layer/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*!
valueB"            2*
(MyModel/Conv1Layer/strided_slice_1/stackГ
*MyModel/Conv1Layer/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*!
valueB"            2,
*MyModel/Conv1Layer/strided_slice_1/stack_1Г
*MyModel/Conv1Layer/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*!
valueB"         2,
*MyModel/Conv1Layer/strided_slice_1/stack_2­
"MyModel/Conv1Layer/strided_slice_1StridedSlicea_in1MyModel/Conv1Layer/strided_slice_1/stack:output:03MyModel/Conv1Layer/strided_slice_1/stack_1:output:03MyModel/Conv1Layer/strided_slice_1/stack_2:output:0*
Index0*
T0*3
_output_shapes!
:         *
ellipsis_mask*
new_axis_mask2$
"MyModel/Conv1Layer/strided_slice_1К
MyModel/Conv1Layer/mulMul#MyModel/Conv1Layer/Reshape:output:0+MyModel/Conv1Layer/strided_slice_1:output:0*
T0*3
_output_shapes!
:         2
MyModel/Conv1Layer/mul═
 MyModel/Conv1Layer/einsum/EinsumEinsumMyModel/Conv1Layer/mul:z:0x_in*
N*
T0*+
_output_shapes
:         *
equationabicf,aif->abc2"
 MyModel/Conv1Layer/einsum/Einsumк
(MyModel/Conv1Layer/MatMul/ReadVariableOpReadVariableOp1mymodel_conv1layer_matmul_readvariableop_resource*
_output_shapes

:*
dtype02*
(MyModel/Conv1Layer/MatMul/ReadVariableOpх
MyModel/Conv1Layer/MatMulBatchMatMulV2x_in0MyModel/Conv1Layer/MatMul/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
MyModel/Conv1Layer/MatMulЙ
MyModel/Conv1Layer/addAddV2)MyModel/Conv1Layer/einsum/Einsum:output:0"MyModel/Conv1Layer/MatMul:output:0*
T0*+
_output_shapes
:         2
MyModel/Conv1Layer/add┼
)MyModel/Conv1Layer/BiasAdd/ReadVariableOpReadVariableOp2mymodel_conv1layer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02+
)MyModel/Conv1Layer/BiasAdd/ReadVariableOp╚
MyModel/Conv1Layer/BiasAddBiasAddMyModel/Conv1Layer/add:z:01MyModel/Conv1Layer/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
MyModel/Conv1Layer/BiasAddњ
MyModel/Conv1Layer/EluElu#MyModel/Conv1Layer/BiasAdd:output:0*
T0*+
_output_shapes
:         2
MyModel/Conv1Layer/Eluм
,MyModel/OutputLayer/Tensordot/ReadVariableOpReadVariableOp5mymodel_outputlayer_tensordot_readvariableop_resource*
_output_shapes

:*
dtype02.
,MyModel/OutputLayer/Tensordot/ReadVariableOpњ
"MyModel/OutputLayer/Tensordot/axesConst*
_output_shapes
:*
dtype0*
valueB:2$
"MyModel/OutputLayer/Tensordot/axesЎ
"MyModel/OutputLayer/Tensordot/freeConst*
_output_shapes
:*
dtype0*
valueB"       2$
"MyModel/OutputLayer/Tensordot/freeъ
#MyModel/OutputLayer/Tensordot/ShapeShape$MyModel/Conv1Layer/Elu:activations:0*
T0*
_output_shapes
:2%
#MyModel/OutputLayer/Tensordot/Shapeю
+MyModel/OutputLayer/Tensordot/GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : 2-
+MyModel/OutputLayer/Tensordot/GatherV2/axisх
&MyModel/OutputLayer/Tensordot/GatherV2GatherV2,MyModel/OutputLayer/Tensordot/Shape:output:0+MyModel/OutputLayer/Tensordot/free:output:04MyModel/OutputLayer/Tensordot/GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2(
&MyModel/OutputLayer/Tensordot/GatherV2а
-MyModel/OutputLayer/Tensordot/GatherV2_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2/
-MyModel/OutputLayer/Tensordot/GatherV2_1/axis╗
(MyModel/OutputLayer/Tensordot/GatherV2_1GatherV2,MyModel/OutputLayer/Tensordot/Shape:output:0+MyModel/OutputLayer/Tensordot/axes:output:06MyModel/OutputLayer/Tensordot/GatherV2_1/axis:output:0*
Taxis0*
Tindices0*
Tparams0*
_output_shapes
:2*
(MyModel/OutputLayer/Tensordot/GatherV2_1ћ
#MyModel/OutputLayer/Tensordot/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2%
#MyModel/OutputLayer/Tensordot/Constл
"MyModel/OutputLayer/Tensordot/ProdProd/MyModel/OutputLayer/Tensordot/GatherV2:output:0,MyModel/OutputLayer/Tensordot/Const:output:0*
T0*
_output_shapes
: 2$
"MyModel/OutputLayer/Tensordot/Prodў
%MyModel/OutputLayer/Tensordot/Const_1Const*
_output_shapes
:*
dtype0*
valueB: 2'
%MyModel/OutputLayer/Tensordot/Const_1п
$MyModel/OutputLayer/Tensordot/Prod_1Prod1MyModel/OutputLayer/Tensordot/GatherV2_1:output:0.MyModel/OutputLayer/Tensordot/Const_1:output:0*
T0*
_output_shapes
: 2&
$MyModel/OutputLayer/Tensordot/Prod_1ў
)MyModel/OutputLayer/Tensordot/concat/axisConst*
_output_shapes
: *
dtype0*
value	B : 2+
)MyModel/OutputLayer/Tensordot/concat/axisћ
$MyModel/OutputLayer/Tensordot/concatConcatV2+MyModel/OutputLayer/Tensordot/free:output:0+MyModel/OutputLayer/Tensordot/axes:output:02MyModel/OutputLayer/Tensordot/concat/axis:output:0*
N*
T0*
_output_shapes
:2&
$MyModel/OutputLayer/Tensordot/concat▄
#MyModel/OutputLayer/Tensordot/stackPack+MyModel/OutputLayer/Tensordot/Prod:output:0-MyModel/OutputLayer/Tensordot/Prod_1:output:0*
N*
T0*
_output_shapes
:2%
#MyModel/OutputLayer/Tensordot/stackЖ
'MyModel/OutputLayer/Tensordot/transpose	Transpose$MyModel/Conv1Layer/Elu:activations:0-MyModel/OutputLayer/Tensordot/concat:output:0*
T0*+
_output_shapes
:         2)
'MyModel/OutputLayer/Tensordot/transpose№
%MyModel/OutputLayer/Tensordot/ReshapeReshape+MyModel/OutputLayer/Tensordot/transpose:y:0,MyModel/OutputLayer/Tensordot/stack:output:0*
T0*0
_output_shapes
:                  2'
%MyModel/OutputLayer/Tensordot/ReshapeЬ
$MyModel/OutputLayer/Tensordot/MatMulMatMul.MyModel/OutputLayer/Tensordot/Reshape:output:04MyModel/OutputLayer/Tensordot/ReadVariableOp:value:0*
T0*'
_output_shapes
:         2&
$MyModel/OutputLayer/Tensordot/MatMulў
%MyModel/OutputLayer/Tensordot/Const_2Const*
_output_shapes
:*
dtype0*
valueB:2'
%MyModel/OutputLayer/Tensordot/Const_2ю
+MyModel/OutputLayer/Tensordot/concat_1/axisConst*
_output_shapes
: *
dtype0*
value	B : 2-
+MyModel/OutputLayer/Tensordot/concat_1/axisА
&MyModel/OutputLayer/Tensordot/concat_1ConcatV2/MyModel/OutputLayer/Tensordot/GatherV2:output:0.MyModel/OutputLayer/Tensordot/Const_2:output:04MyModel/OutputLayer/Tensordot/concat_1/axis:output:0*
N*
T0*
_output_shapes
:2(
&MyModel/OutputLayer/Tensordot/concat_1Я
MyModel/OutputLayer/TensordotReshape.MyModel/OutputLayer/Tensordot/MatMul:product:0/MyModel/OutputLayer/Tensordot/concat_1:output:0*
T0*+
_output_shapes
:         2
MyModel/OutputLayer/Tensordot╚
*MyModel/OutputLayer/BiasAdd/ReadVariableOpReadVariableOp3mymodel_outputlayer_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02,
*MyModel/OutputLayer/BiasAdd/ReadVariableOpО
MyModel/OutputLayer/BiasAddBiasAdd&MyModel/OutputLayer/Tensordot:output:02MyModel/OutputLayer/BiasAdd/ReadVariableOp:value:0*
T0*+
_output_shapes
:         2
MyModel/OutputLayer/BiasAddЋ
MyModel/OutputLayer/EluElu$MyModel/OutputLayer/BiasAdd:output:0*
T0*+
_output_shapes
:         2
MyModel/OutputLayer/Elu}
IdentityIdentity%MyModel/OutputLayer/Elu:activations:0*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         :::::::::::Q M
+
_output_shapes
:         

_user_specified_nameX_in:QM
+
_output_shapes
:         

_user_specified_nameA_in:UQ
/
_output_shapes
:         

_user_specified_nameE_in:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
ќ
Ш
)__inference_Conv1Layer_layer_call_fn_1352
inputs_0
inputs_1
inputs_2
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
identityѕбStatefulPartitionedCall╗
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1inputs_2unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2*
Tout
2*+
_output_shapes
:         **
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8*L
fGRE
C__inference_Conv1Layer_layer_call_and_return_conditional_losses_6332
StatefulPartitionedCallњ
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*|
_input_shapesk
i:         :         :         ::::::::22
StatefulPartitionedCallStatefulPartitionedCall:U Q
+
_output_shapes
:         
"
_user_specified_name
inputs/0:UQ
+
_output_shapes
:         
"
_user_specified_name
inputs/1:YU
/
_output_shapes
:         
"
_user_specified_name
inputs/2:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: 
Ї
Љ
&__inference_MyModel_layer_call_fn_1192
inputs_0
inputs_1
inputs_2
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
identityѕбStatefulPartitionedCallм
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1inputs_2unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8*
Tin
2*
Tout
2*+
_output_shapes
:         *,
_read_only_resource_inputs

	
**
config_proto

CPU

GPU 2J 8*I
fDRB
@__inference_MyModel_layer_call_and_return_conditional_losses_7842
StatefulPartitionedCallњ
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*+
_output_shapes
:         2

Identity"
identityIdentity:output:0*ё
_input_shapess
q:         :         :         ::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:U Q
+
_output_shapes
:         
"
_user_specified_name
inputs/0:UQ
+
_output_shapes
:         
"
_user_specified_name
inputs/1:YU
/
_output_shapes
:         
"
_user_specified_name
inputs/2:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: :

_output_shapes
: :

_output_shapes
: 
ч)
Ў
__inference__traced_save_1451
file_prefix5
1savev2_conv1layer_root_kernel_read_readvariableop.
*savev2_conv1layer_bias_read_readvariableop1
-savev2_outputlayer_kernel_read_readvariableop/
+savev2_outputlayer_bias_read_readvariableop6
2savev2_conv1layer_fgn_0_kernel_read_readvariableop4
0savev2_conv1layer_fgn_0_bias_read_readvariableop6
2savev2_conv1layer_fgn_1_kernel_read_readvariableop4
0savev2_conv1layer_fgn_1_bias_read_readvariableop8
4savev2_conv1layer_fgn_out_kernel_read_readvariableop6
2savev2_conv1layer_fgn_out_bias_read_readvariableop
savev2_1_const

identity_1ѕбMergeV2CheckpointsбSaveV2бSaveV2_1Ј
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*2
StaticRegexFullMatchc
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.part2
ConstЇ
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*<
value3B1 B+_temp_7c7acfaef0c049da87b8563e2c50b890/part2	
Const_1І
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: 2
Selectt

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: 2

StringJoinZ

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :2

num_shards
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 2
ShardedFilename/shardд
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilenameљ
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:
*
dtype0*б
valueўBЋ
B;layer_with_weights-0/root_kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUE2
SaveV2/tensor_namesю
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:
*
dtype0*'
valueB
B B B B B B B B B B 2
SaveV2/shape_and_slicesд
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:01savev2_conv1layer_root_kernel_read_readvariableop*savev2_conv1layer_bias_read_readvariableop-savev2_outputlayer_kernel_read_readvariableop+savev2_outputlayer_bias_read_readvariableop2savev2_conv1layer_fgn_0_kernel_read_readvariableop0savev2_conv1layer_fgn_0_bias_read_readvariableop2savev2_conv1layer_fgn_1_kernel_read_readvariableop0savev2_conv1layer_fgn_1_bias_read_readvariableop4savev2_conv1layer_fgn_out_kernel_read_readvariableop2savev2_conv1layer_fgn_out_bias_read_readvariableop"/device:CPU:0*
_output_shapes
 *
dtypes
2
2
SaveV2Ѓ
ShardedFilename_1/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B :2
ShardedFilename_1/shardг
ShardedFilename_1ShardedFilenameStringJoin:output:0 ShardedFilename_1/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename_1б
SaveV2_1/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*1
value(B&B_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2_1/tensor_namesј
SaveV2_1/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*
valueB
B 2
SaveV2_1/shape_and_slices¤
SaveV2_1SaveV2ShardedFilename_1:filename:0SaveV2_1/tensor_names:output:0"SaveV2_1/shape_and_slices:output:0savev2_1_const^SaveV2"/device:CPU:0*
_output_shapes
 *
dtypes
22

SaveV2_1с
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0ShardedFilename_1:filename:0^SaveV2	^SaveV2_1"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixesг
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix	^SaveV2_1"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

IdentityЂ

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints^SaveV2	^SaveV2_1*
T0*
_output_shapes
: 2

Identity_1"!

identity_1Identity_1:output:0*g
_input_shapesV
T: ::::::::::: 2(
MergeV2CheckpointsMergeV2Checkpoints2
SaveV2SaveV22
SaveV2_1SaveV2_1:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:: 

_output_shapes
::$	 

_output_shapes

:: 


_output_shapes
::

_output_shapes
: "»L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*ф
serving_defaultќ
9
A_in1
serving_default_A_in:0         
=
E_in5
serving_default_E_in:0         
9
X_in1
serving_default_X_in:0         C
OutputLayer4
StatefulPartitionedCall:0         tensorflow/serving/predict:Ўг
л&
layer-0
layer-1
layer-2
layer_with_weights-0
layer-3
layer_with_weights-1
layer-4
	optimizer
loss
_layers
	trainable_variables

	variables
regularization_losses
	keras_api

signatures
*N&call_and_return_all_conditional_losses
O__call__
P_default_save_signature"щ#
_tf_keras_model▀#{"class_name": "Model", "name": "MyModel", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "config": {"name": "MyModel", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 2]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "X_in"}, "name": "X_in", "inbound_nodes": []}, {"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 7]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "A_in"}, "name": "A_in", "inbound_nodes": []}, {"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 7, 2]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "E_in"}, "name": "E_in", "inbound_nodes": []}, {"class_name": "EdgeConditionedConv", "config": {"name": "Conv1Layer", "trainable": true, "dtype": "float32", "channels": 3, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "kernel_constraint": null, "bias_constraint": null, "kernel_network": {"class_name": "__tuple__", "items": [2, 2]}, "root": true}, "name": "Conv1Layer", "inbound_nodes": [[["X_in", 0, 0, {}], ["A_in", 0, 0, {}], ["E_in", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "OutputLayer", "trainable": true, "dtype": "float32", "units": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "OutputLayer", "inbound_nodes": [[["Conv1Layer", 0, 0, {}]]]}], "input_layers": [["X_in", 0, 0], ["A_in", 0, 0], ["E_in", 0, 0]], "output_layers": [["OutputLayer", 0, 0]]}, "build_input_shape": [{"class_name": "TensorShape", "items": [null, 7, 2]}, {"class_name": "TensorShape", "items": [null, 7, 7]}, {"class_name": "TensorShape", "items": [null, 7, 7, 2]}], "is_graph_network": true, "keras_version": "2.3.0-tf", "backend": "tensorflow", "model_config": {"class_name": "Model", "config": {"name": "MyModel", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 2]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "X_in"}, "name": "X_in", "inbound_nodes": []}, {"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 7]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "A_in"}, "name": "A_in", "inbound_nodes": []}, {"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 7, 2]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "E_in"}, "name": "E_in", "inbound_nodes": []}, {"class_name": "EdgeConditionedConv", "config": {"name": "Conv1Layer", "trainable": true, "dtype": "float32", "channels": 3, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "kernel_constraint": null, "bias_constraint": null, "kernel_network": {"class_name": "__tuple__", "items": [2, 2]}, "root": true}, "name": "Conv1Layer", "inbound_nodes": [[["X_in", 0, 0, {}], ["A_in", 0, 0, {}], ["E_in", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "OutputLayer", "trainable": true, "dtype": "float32", "units": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "OutputLayer", "inbound_nodes": [[["Conv1Layer", 0, 0, {}]]]}], "input_layers": [["X_in", 0, 0], ["A_in", 0, 0], ["E_in", 0, 0]], "output_layers": [["OutputLayer", 0, 0]]}}, "training_config": {"loss": null, "metrics": null, "weighted_metrics": null, "loss_weights": null, "sample_weight_mode": null, "optimizer_config": {"class_name": "RMSprop", "config": {"name": "RMSprop", "learning_rate": 0.001, "decay": 0.0, "rho": 0.9, "momentum": 0.0, "epsilon": 1e-07, "centered": false}}}}
ж"Т
_tf_keras_input_layerк{"class_name": "InputLayer", "name": "X_in", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 2]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 2]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "X_in"}}
ж"Т
_tf_keras_input_layerк{"class_name": "InputLayer", "name": "A_in", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 7]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 7]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "A_in"}}
№"В
_tf_keras_input_layer╠{"class_name": "InputLayer", "name": "E_in", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 7, 2]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 7, 7, 2]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "E_in"}}
Ї
kernel_network_layers
root_kernel
bias
trainable_variables
	variables
regularization_losses
	keras_api
*Q&call_and_return_all_conditional_losses
R__call__"╚
_tf_keras_layer«{"class_name": "EdgeConditionedConv", "name": "Conv1Layer", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "Conv1Layer", "trainable": true, "dtype": "float32", "channels": 3, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "kernel_constraint": null, "bias_constraint": null, "kernel_network": {"class_name": "__tuple__", "items": [2, 2]}, "root": true}, "build_input_shape": [{"class_name": "TensorShape", "items": [null, 7, 2]}, {"class_name": "TensorShape", "items": [null, 7, 7]}, {"class_name": "TensorShape", "items": [null, 7, 7, 2]}]}
н

kernel
bias
trainable_variables
	variables
regularization_losses
	keras_api
*S&call_and_return_all_conditional_losses
T__call__"»
_tf_keras_layerЋ{"class_name": "Dense", "name": "OutputLayer", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "OutputLayer", "trainable": true, "dtype": "float32", "units": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 3}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 7, 3]}}
"
	optimizer
 "
trackable_dict_wrapper
C
0
1
2
3
4"
trackable_list_wrapper
f
0
1
2
3
4
5
6
 7
8
9"
trackable_list_wrapper
f
0
1
2
3
4
5
6
 7
8
9"
trackable_list_wrapper
 "
trackable_list_wrapper
╩
!layer_metrics
	trainable_variables
"layer_regularization_losses

	variables

#layers
$non_trainable_variables
%metrics
regularization_losses
O__call__
P_default_save_signature
*N&call_and_return_all_conditional_losses
&N"call_and_return_conditional_losses"
_generic_user_object
,
Userving_default"
signature_map
5
&0
'1
(2"
trackable_list_wrapper
(:&2Conv1Layer/root_kernel
:2Conv1Layer/bias
X
0
1
2
3
4
5
6
 7"
trackable_list_wrapper
X
0
1
2
3
4
5
6
 7"
trackable_list_wrapper
 "
trackable_list_wrapper
Г
)layer_metrics
trainable_variables
*layer_regularization_losses
	variables

+layers
,non_trainable_variables
-metrics
regularization_losses
R__call__
*Q&call_and_return_all_conditional_losses
&Q"call_and_return_conditional_losses"
_generic_user_object
$:"2OutputLayer/kernel
:2OutputLayer/bias
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
Г
.layer_metrics
trainable_variables
/layer_regularization_losses
	variables

0layers
1non_trainable_variables
2metrics
regularization_losses
T__call__
*S&call_and_return_all_conditional_losses
&S"call_and_return_conditional_losses"
_generic_user_object
):'2Conv1Layer/FGN_0/kernel
#:!2Conv1Layer/FGN_0/bias
):'2Conv1Layer/FGN_1/kernel
#:!2Conv1Layer/FGN_1/bias
+:)2Conv1Layer/FGN_out/kernel
%:#2Conv1Layer/FGN_out/bias
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
C
0
1
2
3
4"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
╠

kernel
bias
3trainable_variables
4	variables
5regularization_losses
6	keras_api
*V&call_and_return_all_conditional_losses
W__call__"Д
_tf_keras_layerЇ{"class_name": "Dense", "name": "FGN_0", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "FGN_0", "trainable": true, "dtype": "float32", "units": 2, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 2}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 7, 7, 2]}}
╠

kernel
bias
7trainable_variables
8	variables
9regularization_losses
:	keras_api
*X&call_and_return_all_conditional_losses
Y__call__"Д
_tf_keras_layerЇ{"class_name": "Dense", "name": "FGN_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "FGN_1", "trainable": true, "dtype": "float32", "units": 2, "activation": "relu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 2}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 7, 7, 2]}}
м

kernel
 bias
;trainable_variables
<	variables
=regularization_losses
>	keras_api
*Z&call_and_return_all_conditional_losses
[__call__"Г
_tf_keras_layerЊ{"class_name": "Dense", "name": "FGN_out", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "config": {"name": "FGN_out", "trainable": true, "dtype": "float32", "units": 6, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 2}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 7, 7, 2]}}
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
5
&0
'1
(2"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
Г
?layer_metrics
3trainable_variables
@layer_regularization_losses
4	variables

Alayers
Bnon_trainable_variables
Cmetrics
5regularization_losses
W__call__
*V&call_and_return_all_conditional_losses
&V"call_and_return_conditional_losses"
_generic_user_object
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
Г
Dlayer_metrics
7trainable_variables
Elayer_regularization_losses
8	variables

Flayers
Gnon_trainable_variables
Hmetrics
9regularization_losses
Y__call__
*X&call_and_return_all_conditional_losses
&X"call_and_return_conditional_losses"
_generic_user_object
.
0
 1"
trackable_list_wrapper
.
0
 1"
trackable_list_wrapper
 "
trackable_list_wrapper
Г
Ilayer_metrics
;trainable_variables
Jlayer_regularization_losses
<	variables

Klayers
Lnon_trainable_variables
Mmetrics
=regularization_losses
[__call__
*Z&call_and_return_all_conditional_losses
&Z"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
л2═
@__inference_MyModel_layer_call_and_return_conditional_losses_723
A__inference_MyModel_layer_call_and_return_conditional_losses_1165
@__inference_MyModel_layer_call_and_return_conditional_losses_751
A__inference_MyModel_layer_call_and_return_conditional_losses_1028└
и▓│
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaultsф 
annotationsф *
 
С2р
%__inference_MyModel_layer_call_fn_862
&__inference_MyModel_layer_call_fn_1192
%__inference_MyModel_layer_call_fn_807
&__inference_MyModel_layer_call_fn_1219└
и▓│
FullArgSpec1
args)џ&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaultsџ
p 

 

kwonlyargsџ 
kwonlydefaultsф 
annotationsф *
 
«2Ф
__inference__wrapped_model_517ѕ
І▓Є
FullArgSpec
argsџ 
varargsjargs
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *xбu
sџp
"і
X_in         
"і
A_in         
&і#
E_in         
Ь2в
D__inference_Conv1Layer_layer_call_and_return_conditional_losses_1329б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
М2л
)__inference_Conv1Layer_layer_call_fn_1352б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
№2В
E__inference_OutputLayer_layer_call_and_return_conditional_losses_1383б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
н2Л
*__inference_OutputLayer_layer_call_fn_1392б
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
7B5
!__inference_signature_wrapper_891A_inE_inX_in
е2Цб
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
е2Цб
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
е2Цб
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
е2Цб
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
е2Цб
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 
е2Цб
Ў▓Ћ
FullArgSpec
argsџ
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargsџ 
kwonlydefaults
 
annotationsф *
 њ
D__inference_Conv1Layer_layer_call_and_return_conditional_losses_1329╔ ЉбЇ
ЁбЂ
џ|
&і#
inputs/0         
&і#
inputs/1         
*і'
inputs/2         
ф ")б&
і
0         
џ Ж
)__inference_Conv1Layer_layer_call_fn_1352╝ ЉбЇ
ЁбЂ
џ|
&і#
inputs/0         
&і#
inputs/1         
*і'
inputs/2         
ф "і         Ў
A__inference_MyModel_layer_call_and_return_conditional_losses_1028М
 ЎбЋ
ЇбЅ
џ|
&і#
inputs/0         
&і#
inputs/1         
*і'
inputs/2         
p

 
ф ")б&
і
0         
џ Ў
A__inference_MyModel_layer_call_and_return_conditional_losses_1165М
 ЎбЋ
ЇбЅ
џ|
&і#
inputs/0         
&і#
inputs/1         
*і'
inputs/2         
p 

 
ф ")б&
і
0         
џ І
@__inference_MyModel_layer_call_and_return_conditional_losses_723к
 їбѕ
ђб}
sџp
"і
X_in         
"і
A_in         
&і#
E_in         
p

 
ф ")б&
і
0         
џ І
@__inference_MyModel_layer_call_and_return_conditional_losses_751к
 їбѕ
ђб}
sџp
"і
X_in         
"і
A_in         
&і#
E_in         
p 

 
ф ")б&
і
0         
џ ы
&__inference_MyModel_layer_call_fn_1192к
 ЎбЋ
ЇбЅ
џ|
&і#
inputs/0         
&і#
inputs/1         
*і'
inputs/2         
p

 
ф "і         ы
&__inference_MyModel_layer_call_fn_1219к
 ЎбЋ
ЇбЅ
џ|
&і#
inputs/0         
&і#
inputs/1         
*і'
inputs/2         
p 

 
ф "і         с
%__inference_MyModel_layer_call_fn_807╣
 їбѕ
ђб}
sџp
"і
X_in         
"і
A_in         
&і#
E_in         
p

 
ф "і         с
%__inference_MyModel_layer_call_fn_862╣
 їбѕ
ђб}
sџp
"і
X_in         
"і
A_in         
&і#
E_in         
p 

 
ф "і         Г
E__inference_OutputLayer_layer_call_and_return_conditional_losses_1383d3б0
)б&
$і!
inputs         
ф ")б&
і
0         
џ Ё
*__inference_OutputLayer_layer_call_fn_1392W3б0
)б&
$і!
inputs         
ф "і         з
__inference__wrapped_model_517л
 ѓб
xбu
sџp
"і
X_in         
"і
A_in         
&і#
E_in         
ф "=ф:
8
OutputLayer)і&
OutputLayer         ї
!__inference_signature_wrapper_891Т
 ўбћ
б 
їфѕ
*
A_in"і
A_in         
.
E_in&і#
E_in         
*
X_in"і
X_in         "=ф:
8
OutputLayer)і&
OutputLayer         