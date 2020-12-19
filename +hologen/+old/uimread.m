function [I,map]=uimread()
%%读入图片函数
I=0;map=0;
global old_pr
if old_pr==0
    old_pr=[];
end
[fn,pn]=uigetfile({'*.bmp;*.jpg; *.tiff;*.tif; *.gif; *.png','图像文件';'*.*','所有文件'},'打开图片',old_pr,'MultiSelect','on');
old_pr=pn;
if isequal(fn,0)||isequal(pn,0)
    return;
end
if iscell(fn)
    [~,N]=size(fn);
    [Is,map]=imread([pn fn{1}]);
        [m,n,~]=size(Is);
        I=zeros(m,n,N);
        I=uint8(I);
        I(:,:,1)=Is(:,:,1);
    for i=2:N
        [Is,map]=imread([pn fn{i}]);
        I(:,:,i)=Is(:,:,1);
    end
else
[Is,map]=imread([pn fn]);
I=Is(:,:,1);
end