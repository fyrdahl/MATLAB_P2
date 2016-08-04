
relaxParams(1) = 220;
relaxParams(2) = 200;
relaxParams(3) = 300;
relaxParams(4) = 280;
x = fminsearch(@(params) dot_product_sim_SE(params,relaxParams),[0,0,90,100])
%%
[x,fval,exitflag,output]=fminsearchbnd(@(params) dot_product_sim_SE(params,relaxParams),[10,100,90,100],[0 0 60 150],[1000 1000 120 210])