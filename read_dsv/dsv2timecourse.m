function output_data = dsv2timecourse ( dsv )
% Take a Siemens dsv as input
%
% output a timecourse for the file (basically unpack the DSV)
%
% MDR April 2008
%
%

input_data = dsv.output_data;

look_for_repeats = 1;
output_data(1) = input_data(1);
out_counter = 2;
for counter = 2 : length(input_data)
    
    if (look_for_repeats == -1)
        look_for_repeats = look_for_repeats + 1;
    elseif (look_for_repeats == 0)
        output_data(out_counter) = output_data(out_counter-1) + input_data(counter);
        out_counter = out_counter + 1;
        look_for_repeats = look_for_repeats + 1;
    else
        output_data(out_counter) = output_data(out_counter-1) + input_data(counter);
        out_counter = out_counter + 1;
        if (input_data(counter) == input_data(counter-1))
            for extra_points = 1 : input_data(counter+1)
                output_data(out_counter) = output_data(out_counter-1) + input_data(counter);
                out_counter = out_counter + 1;
            end
            look_for_repeats = -1; % Jump over this value
        end
    end
    
end


