function gammahosei=gammahosei(g_norm,height,width)

gammahosei=zeros(height*width,3);
for n=1:height*width
    for k=1:3
        if g_norm(n,k)>0.00313
            gammahosei(n,k)=1.055.*g_norm(n,k)^(1/2.4)-0.055;
        else
            gammahosei(n,k)=12.92.*g_norm(n,k);
        end
    end
end

end