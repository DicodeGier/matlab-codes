clear
clc
%%input 
input_matrix = [7,2,8,3;3,8,4,6;2,3,7,7];
sizes = size(input_matrix);
rows = sizes(1);
columns = sizes(2);

%%initialize the output
names = ["xxxx",'1xxx','2xxx','3xxx','4xxx','12xx','13xx','14xx','21xx','23xx','24xx','31xx','32xx','34xx','41xx','42xx','43xx','1234','1243','1324','1342','1423','1432','2134','2143','2314','2341','2413','2431','3124','3142','3214','3241','3412','3421','4123','4132','4213','4231','4312','4321'];
output = cell(length(names),3);
for i = 1:length(names)
    output{i,1} = names(i);
end
%%%%%first iteration: no constraints on jobs
%%calculate rows
first_row = sum(input_matrix(1,:)) + min(sum(input_matrix(2:end,:)));
second_row = min(input_matrix(1,:)) + sum(input_matrix(2,:)) + min(input_matrix(3,:));
third_row = min(sum(input_matrix(1:2,:))) + sum(input_matrix(3,:));
%%put rows in struct and determine max row
output_struct = struct();
output_struct.xxxx.first_row = first_row;
output_struct.xxxx.second_row = second_row;
output_struct.xxxx.third_row = third_row;
output_value = max([first_row, second_row, third_row]);
%%fill in cell
output{1,1} = 'xxxx';
output{1,2} = output_struct.xxxx;
output{1,3} = output_value;

%%%%%now determine ixxx with i for every column/job
for j = 1:columns
    first_row = sum(input_matrix(1,:)) + min(sum(input_matrix(2:end,[1:j-1 j+1:end])));
    second_row = input_matrix(1,j) + sum(input_matrix(2,:)) + min(input_matrix(3,[1:j-1 j+1:end]));
    third_row = sum(input_matrix(:,j)) + sum(input_matrix(3,[1:j-1 j+1:end]));
    output_struct.jxxx.first_row = first_row;
    output_struct.jxxx.second_row = second_row;
    output_struct.jxxx.third_row = third_row;
    output{j+1,2} = output_struct.jxxx;
    output_value = max([first_row, second_row, third_row]);
    output{j+1,3} = output_value;
end

%%%%%now determine all combinations in which two jobs are scheduled
k = 6;
for i = 1:columns
    for j = find([1:columns]~=i)
        submatrix = input_matrix(:,[1:4]~=i);
        if i < j
            submatrix = submatrix(:,[1:3]~=j-1);
        else
            submatrix = input_matrix(:,[1:3]~=j);
        end
        first_row = sum(input_matrix(1,:)) + min(sum(submatrix(2:end,:)));
        second_row = max([input_matrix(1,i) + input_matrix(1,j), input_matrix(1,i) + input_matrix(2,i)]) + input_matrix(2,j) + sum(submatrix(2,:)) + min(submatrix(3,:));
        third_row = max([sum(input_matrix(:,i)), second_row - sum(submatrix(2,:)) - min(submatrix(3,:))]) + input_matrix(3,j) + sum(submatrix(3,:));
        output_struct.ijxx.first_row = first_row;
        output_struct.ijxx.second_row = second_row;
        output_struct.ijxx.third_row = third_row;
        output{k,2} = output_struct.ijxx;
        output_value = max([first_row, second_row, third_row]);
        output{k,3} = output_value; 
        k = k+1;
    end
end

%%%%% now determine all combinations in which the job order is fully
%%%%% specified
z = 18;
for i = 1:4
    for j = 1:4
        if j == i
            continue
        else
        for k = 1:4
            if k == j || k == i
                continue
            else
            for l = 1:4
                if l == j || l == k || l == i
                    continue
                else
                     ji_1 = input_matrix(1,i);
                     ji_2 = ji_1 + input_matrix(2,i);
                     ji_3 = ji_2 + input_matrix(3,i);
                     
                     jj_1 = ji_1 + input_matrix(1,j);
                     jj_2 = max([jj_1, ji_2]) + input_matrix(2,j);
                     jj_3 = max([jj_2, ji_3]) + input_matrix(3,j);
                     
                     jk_1 = jj_1 + input_matrix(1,k);
                     jk_2 = max([jk_1, jj_2]) + input_matrix(2,k);
                     jk_3 = max([jk_2, jj_3]) + input_matrix(3,k);
                     
                     jl_1 = jk_1 + input_matrix(1,l);
                     jl_2 = max([jl_1, jk_2]) + input_matrix(2,l);
                     jl_3 = max([jl_2, jk_3]) + input_matrix(3,l);
                     
                     first_row = jl_1;
                     second_row = jl_2;
                     third_row = jl_3;
                     
                     output_struct.ijkl.first_row = first_row;
                     output_struct.ijkl.second_row = second_row;
                     output_struct.ijkl.third_row = third_row;
                     
                     output{z,2} = output_struct.ijkl;
                     output_value = max([first_row, second_row, third_row]);
                     output{z,3} = output_value;
                     z = z+1;
                end
            end
            end
        end
        end
    end
end

output_subset = output(18:end,:);
output_mat = cell2mat(output_subset(:,3));
min_solution = min(output_mat)
row_index = find(output_mat(:,1) == min_solution);
min_job_order = output_subset{row_index,1}


