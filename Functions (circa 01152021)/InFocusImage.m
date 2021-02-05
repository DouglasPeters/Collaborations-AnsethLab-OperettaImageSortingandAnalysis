

function [Ch1,Ch2,Ch3,Ch4,Focus_Ch1,Focus_Ch2,Focus_Ch3,Focus_Ch4,ZFocus] = InFocusImage(I,ZPlanes,ResY,ResX,Channels,f)

for m = 1:ZPlanes
    Ch1(:,:,m) = I{1,1}{m,1}(:,:);
   
    Ch1_grad(m).GradImage(:,:) = imgradient(Ch1(:,:,m));
    Ch1_grad(m).Mean = mean(Ch1_grad(m).GradImage(:,:),'all');
    Ch1_grad(m).Percentile = prctile(Ch1_grad(m).GradImage(:,:),[95 99],'all');
    Filtered(:,:,m) = Ch1_grad(m).GradImage(Ch1_grad(m).Percentile(1,1) < Ch1_grad(m).GradImage < Ch1_grad(m).Percentile(2,1));    
    FilteredMean(m,1) = mean(mean(Filtered(:,:,m)));
    
    Ch2(:,:,m) = I{1,1}{ZPlanes+m,1}(:,:);
    if Channels>2, Ch3(:,:,m) = I{1,1}{(2*ZPlanes)+m,1}(:,:); 
    else Ch3 = zeros(1); end;
    if Channels>3, Ch4(:,:,m) = I{1,1}{(3*ZPlanes)+m,1}(:,:); 
    else Ch4 = zeros(1); end;
end
[MAX(f,1),ZFocus(f,1)] = max(FilteredMean);

if ZFocus(f,1) == 1 && ZPlanes > 1
    for y = 1:ResY
        for x = 1:ResX
            Focus_Ch1(y,x) = max([Ch1(y,x,ZFocus(f,1)),Ch1(y,x,(ZFocus(f,1)+1))]);
            Focus_Ch2(y,x) = max([Ch2(y,x,ZFocus(f,1)),Ch2(y,x,(ZFocus(f,1)+1))]);
            if Channels>2, Focus_Ch3(y,x) = max([Ch3(y,x,ZFocus(f,1)),Ch3(y,x,(ZFocus(f,1)+1))]); else end;
            if Channels>3, Focus_Ch4(y,x) = max([Ch4(y,x,ZFocus(f,1)),Ch4(y,x,(ZFocus(f,1)+1))]); else end;
        end
    end
elseif ZFocus(f,1) == ZPlanes && ZPlanes > 1
    for y = 1:ResY
        for x = 1:ResX
            Focus_Ch1(y,x) = max([Ch1(y,x,ZFocus(f,1)),Ch1(y,x,(ZFocus(f,1)-1))]);
            Focus_Ch2(y,x) = max([Ch2(y,x,ZFocus(f,1)),Ch2(y,x,(ZFocus(f,1)-1))]);
            if Channels>2, Focus_Ch3(y,x) = max([Ch3(y,x,ZFocus(f,1)),Ch3(y,x,(ZFocus(f,1)-1))]); else end;
            if Channels>3, Focus_Ch4(y,x) = max([Ch4(y,x,ZFocus(f,1)),Ch4(y,x,(ZFocus(f,1)-1))]); else end;
        end
    end
elseif ZPlanes > 1
    for y = 1:ResY
        for x = 1:ResX
            Focus_Ch1(y,x) = max([Ch1(y,x,ZFocus(f,1)),Ch1(y,x,(ZFocus(f,1)+1)),Ch1(y,x,(ZFocus(f,1)-1))]);
            Focus_Ch2(y,x) = max([Ch2(y,x,ZFocus(f,1)),Ch2(y,x,(ZFocus(f,1)+1)),Ch2(y,x,(ZFocus(f,1)-1))]);
            if Channels>2, Focus_Ch3(y,x) = max([Ch3(y,x,ZFocus(f,1)),Ch3(y,x,(ZFocus(f,1)+1)),Ch3(y,x,(ZFocus(f,1)-1))]); else end;
            if Channels>3, Focus_Ch4(y,x) = max([Ch4(y,x,ZFocus(f,1)),Ch4(y,x,(ZFocus(f,1)+1)),Ch4(y,x,(ZFocus(f,1)-1))]); else end;
        end
    end
else
    Focus_Ch1 = Ch1;
    Focus_Ch2 = Ch2;
    if Channels>2, Focus_Ch3 = Ch3; else end
    if Channels>3, Focus_Ch4 = Ch4; else end
end

if exist('Focus_Ch3') == 0, Focus_Ch3 = Ch3; else end;
if exist('Focus_Ch4') == 0, Focus_Ch4 = Ch4; else end;

end