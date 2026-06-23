function dy = net_int(omega,bias,win,sig,time,t,y)
if max(abs(sig))>0
u  = interp1(time,sig,t,'nearest');
else 
u = 0;
end

dy = -y + tanh(omega*y+bias+win*u);
end