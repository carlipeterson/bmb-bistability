function F = mySystemTick(vars,params)

    F = vars - tickMap(vars,params);

end