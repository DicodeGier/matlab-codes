clear

load('K_paths_values_10_paths.mat')
cities = cell(49,1);
cities{1} = 'Amsterdam';
cities{2} = 'Berlijn';
cities{3} = 'Hamburg';
cities{4} = 'München';
cities{5} = 'Frankfurt';
cities{6} = 'Hannover';
cities{7} = 'Keulen';
cities{8} = 'Brussel';
cities{9} = 'Lille';
cities{10} = 'Kopenhagen';
cities{11} = 'Oslo';
cities{12} = 'Stockholm';
cities{13} = 'Parijs';
cities{14} = 'Lyon';
cities{15} = 'Marseille';
cities{16} = 'Straatsburg';
cities{17} = 'Bordeaux';
cities{18} = 'Perpignan';
cities{19} = 'Bern';
cities{20} = 'Barcelona';
cities{21} = 'Madrid';
cities{22} = 'Sevilla';
cities{23} = 'Burgos';
cities{24} = 'Cordoba';
cities{25} = 'Praag';
cities{26} = 'Milaan';
cities{27} = 'Rome';  
cities{28} = 'Napels';
cities{29} =    'Venetië';
cities{30} =    'Bologna';
cities{31} =    'Florence';
cities{32} =    'Linz';
cities{33} =    'Wenen';
cities{34} =    'Ljubljana';
cities{35} =    'Zagreb'; 
cities{36} =    'Boedapest'; 
cities{37} =    'Warschau'; 
cities{38} =    'Katowice';
cities{39} =    'Belgrado';
cities{40} =    'Podgorica';
cities{41} =    'Boekarest';
cities{42} =    'Sofia';
cities{43} =    'Skopje';
cities{44} =    'Thessaloniki';
 cities{45} =   'Athene';
 cities{46} =   'Istanbul';
 cities{47} =   'Londen';
 cities{48} =   'Birmingham';
 cities{49} =   'Glasgow';
 
 total_edges_shsr = sum(sum(total_upgrades_shsr_10 > 0));
 cell_edges_shsr = cell(total_edges_shsr,3);
 counter = 1;
 for i = 1:49
     for j = 1:49
     if total_upgrades_shsr_10(i,j) > 0
        cell_edges_shsr{counter,1} = cities{i};
        cell_edges_shsr{counter,2} = cities{j};
        cell_edges_shsr{counter,3} = total_upgrades_shsr_10(i,j);
        counter = counter + 1;
     end
     end
 end
 
 sorted_cell_edges_shsr = sortrows(cell_edges_shsr,3);
 
 total_edges_hsr = sum(sum(total_upgrades_hsr_10 > 0));
 cell_edges_hsr = cell(total_edges_hsr,3);
 counter = 1;
 for i = 1:49
     for j = 1:49
     if total_upgrades_hsr_10(i,j) > 0
        cell_edges_hsr{counter,1} = cities{i};
        cell_edges_hsr{counter,2} = cities{j};
        cell_edges_hsr{counter,3} = total_upgrades_hsr_10(i,j);
        counter = counter + 1;
     end
     end
 end
 
 sorted_cell_edges_hsr = sortrows(cell_edges_hsr,3);