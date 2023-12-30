function Y=voice(x,f) %更改采样率使基频改变 f>1降低（女变男）;f<1升高（男变女）
    f=round(f*1000);    
    d=resample(x,f,1000); %时长整合使语音文件恢复原来时长    
    W=400;  
    Wov=W/2;    
    Kmax=W*2;   
    Wsim=Wov;    
    xdecim=8;    
    kdecim=2; X=d';    
    F=f/1000;   
    Ss =W-Wov;    
    xpts = size(X,2);    
    ypts = round(xpts / F);    
    Y = zeros(1, ypts);    
    xfwin = (1:Wov)/(Wov+1);    
    ovix = (1-Wov):0; newix = 1:(W-Wov);    
    simix = (1:xdecim:Wsim) - Wsim;    
    padX = [zeros(1, Wsim), X, zeros(1,Kmax+W-Wov)];   
    Y(1:Wsim) = X(1:Wsim); lastxpos = 0; km = 0;   
    for ypos = Wsim:Ss:(ypts-W)        
        xpos = round(F * ypos);        
        kmpred = km + (xpos - lastxpos);        
        lastxpos = xpos;        
        if (kmpred <= Kmax)        
            km = kmpred;       
        else        
            ysim = Y(ypos + simix);        
            rxy = zeros(1, Kmax+1);        
            rxx = zeros(1, Kmax+1);        
            Kmin = 0;        
            for k = Kmin:kdecim:Kmax        
                xsim = padX(Wsim + xpos + k + simix);            
                rxx(k+1) = norm(xsim);            
                rxy(k+1) = (ysim * xsim');        
            end        
            Rxy = (rxx ~= 0).*rxy./(rxx+(rxx==0));        
            km = min(find(Rxy == max(Rxy))-1);        
        end        
        xabs = xpos+km;       
        Y(ypos+ovix) = ((1-xfwin).*Y(ypos+ovix)) + (xfwin.*padX(Wsim+xabs+ovix));        
        Y(ypos+newix) = padX(Wsim+xabs+newix);    
    end
end
