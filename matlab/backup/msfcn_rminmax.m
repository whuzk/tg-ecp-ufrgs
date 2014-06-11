function msfcn_rminmax(block)
setup(block);

function setup(block)
% Register number of dialog parameters 
block.NumDialogPrms = 2;
% Register number of ports
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;
% Override input port properties
block.InputPort(1).Dimensions = 1;
block.InputPort(1).DatatypeID = 6; % int32
block.InputPort(1).Complexity = 'Real';
block.InputPort(1).SamplingMode = 'Sample';
block.InputPort(1).DirectFeedthrough = true;
% Override output port properties
block.OutputPort(1).Dimensions = 1;
block.OutputPort(1).DatatypeID = 6; % int32
block.OutputPort(1).Complexity = 'Real';
block.OutputPort(1).SamplingMode = 'Sample';
% Register sample times (Inherited)
block.SampleTimes = [-1 0];
% Specify the block simStateCompliance
block.SimStateCompliance = 'DefaultSimState';
% Register methods
block.RegBlockMethod('CheckParameters', @CheckParam);
block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs); % Required
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Terminate', @Terminate); % Required

function CheckParam(block)
wsize = block.DialogPrm(1).Data;
if wsize <= 0
   error('The window size must be greater than zero.');
end

function DoPostPropSetup(block)
wsize = block.DialogPrm(1).Data;
block.NumDworks = 5;
% State variable #1
block.Dwork(1).Name = 'val';
block.Dwork(1).Dimensions = wsize;
block.Dwork(1).DatatypeID = 6; % int32
block.Dwork(1).Complexity = 'Real';
block.Dwork(1).UsedAsDiscState = true;
% State variable #2
block.Dwork(2).Name = 'pos';
block.Dwork(2).Dimensions = wsize;
block.Dwork(2).DatatypeID = 6; % int32
block.Dwork(2).Complexity = 'Real';
block.Dwork(2).UsedAsDiscState = true;
% State variable #3
block.Dwork(3).Name = 'first';
block.Dwork(3).Dimensions = 1;
block.Dwork(3).DatatypeID = 6; % int32
block.Dwork(3).Complexity = 'Real';
block.Dwork(3).UsedAsDiscState = true;
% State variable #4
block.Dwork(4).Name = 'count';
block.Dwork(4).Dimensions = 1;
block.Dwork(4).DatatypeID = 6; % int32
block.Dwork(4).Complexity = 'Real';
block.Dwork(4).UsedAsDiscState = true;
% State variable #5
block.Dwork(5).Name = 'index';
block.Dwork(5).Dimensions = 1;
block.Dwork(5).DatatypeID = 6; % int32
block.Dwork(5).Complexity = 'Real';
block.Dwork(5).UsedAsDiscState = true;

function Start(block)
wsize = block.DialogPrm(1).Data;
block.Dwork(1).Data = zeros(1,wsize,'int32');
block.Dwork(2).Data = zeros(1,wsize,'int32');
block.Dwork(3).Data = int32(1);
block.Dwork(4).Data = int32(0);
block.Dwork(5).Data = int32(1);

function Outputs(block)
block.OutputPort(1).Data = block.Dwork(1).Data(block.Dwork(3).Data);

function Update(block)
% Get input, state variables and block parameters
x = block.InputPort(1).Data;
first = block.Dwork(3).Data;
count = block.Dwork(4).Data;
i = block.Dwork(5).Data;
wsize = block.DialogPrm(1).Data;
ismax = block.DialogPrm(2).Data;
% Main routine
j = count;
if ismax
    while j > 0 && x >= block.Dwork(1).Data(mod(first+j-2,wsize)+1)
        j = j - 1;
    end
else
    while j > 0 && x <= block.Dwork(1).Data(mod(first+j-2,wsize)+1)
        j = j - 1;
    end
end
idx = mod(first+j-1,wsize)+1;
count = j + 1;
if count > wsize || block.Dwork(2).Data(first) <= i-wsize
    count = count - 1;
    first = mod(first,wsize)+1;
end
% Update state variables
block.Dwork(1).Data(idx) = x;
block.Dwork(2).Data(idx) = i;
block.Dwork(3).Data = first;
block.Dwork(4).Data = count;
block.Dwork(5).Data = i + 1;

function Terminate(block)