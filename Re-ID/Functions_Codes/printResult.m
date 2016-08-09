function printResult(OC, JPM, data,time,rank)
 

disp(' DATA | Method |     Recongnition Rate     | time')
disp(['      |        | Rank: ',num2str(rank), ' | sec'])
disp([data,' | Org.CO |',num2str(round(OC(rank)'*10)/10),' | 0.00'])
disp([data,' | ourJPM | ',num2str(round(JPM(rank)'*10)/10), ' | ', num2str(round(10*time)/10)])
