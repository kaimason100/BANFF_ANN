function dy = net_int(eta,phi,bias,win,sig,time,t,y)
if max(abs(sig))>0
u  = interp1(time,sig,t,'nearest');
else 
u = 0;
end

dy = -y + phi'*(tanh(eta*y+bias+win*u));
end