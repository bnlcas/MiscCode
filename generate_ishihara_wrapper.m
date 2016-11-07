function [] = generate_ishihra_plots_wrapper()
%% This function generates a series of ishihara prints

root_dir = '/Users/changlab/Documents/changrepo/matlab/analysis/Ishihara';
save_folder = [root_dir '/ishihara_images'];
templates_folder = [root_dir '/templates'];

if exist(save_folder) ~= 7
    mkdir(save_folder)
end


%% get template files
template_file_type = 'png';
templates = dir([templates_folder, '/*.' template_file_type]);


%% Loop Through, Randomly Selecting Templates and color patterns:
iterations = 20;
for i = 1:iterations
    fig = figure;
    hold on;
    % Randomly select a template:
    template_name = templates(randperm(length(templates),1)).name;

    % Randomly Select RGB Triplets:
    %tmp = (0.7*rand+(1-0.7)/2);
    %Figure_RGB = [rand rand rand];
    %Ground_RGB = [rand rand rand];
    Figure_RGB = [(0.8*rand+(1-0.8)/2) (0.8*rand+(1-0.8)/2) (0.8*rand+(1-0.8)/2)];
    Ground_RGB = abs(Figure_RGB + [(0.5*rand-0.25) (0.5*rand-0.25) (0.5*rand-0.25)]);
    
    var = [0.07 0.07 0.07];


    generate_ishihara([templates_folder filesep template_name], Figure_RGB, Ground_RGB, var);
    axis off;
    
    ishihara_name = strrep([template_name(1:(end-4)) ' FIG ' num2str(round(Figure_RGB*100)) ' GND ' num2str(round(Ground_RGB*100))],' ', '_');

    print(fig, [save_folder '/' ishihara_name], '-dpng', '-r300'); % may be better as a .svg file...
    close(fig);
    
end



end

