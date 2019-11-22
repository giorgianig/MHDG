function bool = structCompare(s1,s2)

global messages
messages = {};
global rootStructNames
rootStructNames = {};

% if isempty(messages)
bool = structCompareImpl(s1,s2);
if bool
    disp('compared structs are equal.')
else
    disp(sprintf('compared structs are not equal.\n'))
    n=numel(messages);
    for i=n:-1:1
        disp(messages{i});
    end
end

end

function bool=structCompareImpl(s1,s2)

global rootStructNames

names1=fieldnames(s1);
names2=fieldnames(s2);

size1=numel(names1);
size2=numel(names2);

bool = true;

if size1==0 || size2==0
    error('empty structs.');
    bool = false;    
else
    for i=size1:-1:1
       field1=names1{i};
       field2=[];
       for j=1:size2
           if(strcmp(field1,names2{j}))
               field2=names2{j};
               break;
           end
       end
       if ~strcmp(field1,field2)
           disp_message_FieldNotFound(field1,'arg1','arg2');
           bool = false;
           s1=rmfield(s1,field1);           
       else
           content1=getfield(s1,field1);
           content2=getfield(s2,field2);
           classContent1=class(content1);
           classContent2=class(content2);
           
           if ~strcmp(classContent1,classContent2)
                disp_message_DifferentClasses(field1);
                bool = false;                
           else
               if strcmp(classContent1,'struct')
                   rootStructNames = [rootStructNames field1];
                   if ~structCompareImpl(content1,content2)
                         disp_message_DifferentFieldContents(field1,'struct');
                         bool = false;
                   end
                   
               elseif strcmp(classContent1,'char')
                   if ~strcmp(content1,content2)
                       disp_message_DifferentFieldContents(field1,'char');
                       bool = false;
                   end
               elseif strcmp(classContent1,'double')
                   if isscalar(content1)
                       if content1~=content2
                           disp_message_DifferentFieldContents(field1,'double scalar');
                           bool = false;                           
                       end
                   else
                       if ~isequal(content1,content2)
                           disp_message_DifferentFieldContents(field1,'double array');                      
                           bool = false; 
                       end
                   end
               elseif strcmp(classContent1,'int32')
                   if isscalar(content1)
                       if content1~=content2
                           disp_message_DifferentFieldContents(field1,'int32 scalar');                       
                           bool = false;                           
                       end
                   else
                       if ~isequal(content1,content2)
                           disp_message_DifferentFieldContents(field1,'int32 array');
                           bool = false;                           
                       end
                   end
               else
%                   bool = false;
%                   error('class type %s of field %s cannot be analyzed.',classContent1,field1);
					disp_message_CellsWarning(field1);
               end               
               
           end
           
           s1=rmfield(s1,field1);
           s2=rmfield(s2,field2);
       end
   
    end
    
    remainedNames2=fieldnames(s2);
    if numel(remainedNames2)>0
         bool = false;
         for i=1:numel(remainedNames2)
             disp_message_FieldNotFound(remainedNames2{i},'arg2','arg1');
         end
    end

    if size1~=size2
    disp_message_0();
    bool = false;
    end
    
    if numel(rootStructNames)>0
        rootStructNames=rootStructNames(1:(numel(rootStructNames)-1));
    end

end
end

function [] = disp_message_0()

[root,spaceTab] = getRoot();
rootNotDot = root(1:end-1);

if strcmp(root,'')
    text=sprintf(['different number of fields in structs.']);
else
    text=sprintf(['different number of fields in structs [' rootNotDot ']']);
end

addMessage([spaceTab text]);

end

function [] = disp_message_FieldNotFound(aField,s1Name,s2Name)

[root,spaceTab] = getRoot();

fullField=[root aField];
text=sprintf('field [%s] of %s not found in %s',fullField,s1Name,s2Name);

addMessage([spaceTab text]);

end

function [] = disp_message_DifferentClasses(aField)

[root,spaceTab] = getRoot();

fullField=[root aField];
text=sprintf('different class fields [%s]',fullField);

addMessage([spaceTab text]);

end
                      
function [] = disp_message_DifferentFieldContents(aField,aClass)

[root,spaceTab] = getRoot();

fullField=[root aField];
text=sprintf('contents of field [%s],typed as -%s-, are not matching',fullField,aClass);

addMessage([spaceTab text]);

end

function [] = disp_message_CellsWarning(aField)

[root,spaceTab] = getRoot();

fullField=[root aField];
text=sprintf('Warning: imposible to compare cells. Field [%s],typed as -cell-',fullField);

addMessage([spaceTab text]);

end


function [root,spaceTab] = getRoot()

global rootStructNames
root='';
if numel(rootStructNames)>0
    root=rootStructNames{1};
    for i=2:numel(rootStructNames)
        root=[root '.' rootStructNames{i}];
    end
    root=[root '.'];
end

spaceTab='';
if numel(rootStructNames)>0
    for i=1:numel(rootStructNames)
        spaceTab=[spaceTab '    '];
    end
end 

end

function []=addMessage(text)

global messages
messages=[messages text];

end
