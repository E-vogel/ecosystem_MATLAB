clear
close all

% Initial population of each agent
n_grass = 200;
n_sheep = 100;
n_dog = 15;

% Step size and time
h = 0.01;
t = 0:h:10;

% Initial state (position and velocity)
rang = 10; 

X_grass = rang*(2*rand(n_grass,2) - 1); % grass initial position
X_sheep = 2*(2*rand(n_sheep,2) - 1); % sheep initial position
X_dog = rang*rand(n_dog,2) - rang; % dog initial position

V_sheep = 2*rand(n_sheep,2) - 1; % sheep initial velocity
V_dog = 2*rand(n_dog,2) - 1; % dog initial velocity

Vs = 3*(2*rand(1,2) - 1); % sheep moving speed

% Lifespan settings
lifespanmax = 300; % maximum lifespan
lifespan_sheep = lifespanmax*ones(n_sheep,1);
lifespan_dog = randi([50 300],n_dog,1);

% Hunger settings
hungermax = 20; % state of satiety
hunger_sheep = hungermax*ones(n_sheep,1);
hunger_dog = randi([1 10],n_dog,1);

% Gains
k_sheep_av = 0.01;
k_sheep_rp = 1.2;
k_sheep_at = 0.2;

k_dog_at = 0.15;

% Storage for values
n_sheep_str = zeros(length(t),1);
n_dog_str = zeros(length(t),1);

n_sheep_str(1) = n_sheep;
n_dog_str(1) = n_dog;

X_grass_str = cell(length(t),1); 
X_sheep_str = cell(length(t),1); 
X_dog_str = cell(length(t),1); 

X_grass_str{1} = X_grass;
X_sheep_str{1} = X_sheep;
X_dog_str{1} = X_dog;

% Simulation
disp("Give me a second...")
for step = 2:length(t)
    disp([step n_sheep n_dog])
    % Update velocity
    V_sheep = V_sheep - k_sheep_av*(V_sheep - sum(V_sheep)/n_sheep); % Void (alignment)
    V_sheep = V_sheep - k_sheep_av*(V_sheep - sum(X_sheep)/n_sheep); % Void (cohesion)

    if mod(step,60) == 0 % Direction of travel (cycle 60)
        Vs = 3*(2*rand(1,2) - 1);
    end

    V_sheep = V_sheep + 3*(2*rand(size(V_sheep)) - 1) + ones(size(V_sheep)).*Vs; % Random + direction
    V_dog = V_dog + 10*(2*rand(size(V_dog)) - 1); % Random

    for i = 1:n_sheep
        for j = 1:n_dog
            distance_tmp = norm(X_sheep(i,:)-X_dog(j,:));
            if distance_tmp < rang && hunger_dog(j) <= 0 % Starvation
                V_dog(j,:) = V_dog(j,:) - k_dog_at*(X_dog(j,:) - X_sheep(i,:))/distance_tmp;
            end
            if distance_tmp < rang*0.4 % Avoidance
                V_sheep(i,:) = V_sheep(i,:) + k_sheep_rp*(X_sheep(i,:) - X_dog(j,:))/distance_tmp;
            end
        end
        for ii = 1:n_sheep
            distance_tmp = norm(X_sheep(i,:)-X_sheep(ii,:));
            if i ~= ii && distance_tmp > rang*0.6 && distance_tmp < rang*1.2 % Approach the same species
                V_sheep(i,:) = V_sheep(i,:) - k_sheep_at*(X_sheep(i,:) - X_sheep(ii,:))/distance_tmp;
            end
        end
        
    end
    V_sheep = V_sheep*0.9; % Deceleration
    V_dog = V_dog*0.9;

    % Wall boundary conditions
%     V_dog = V_dog.*(2*(abs(X_dog) < rang) - 1);
%     V_sheep = V_sheep.*(2*(abs(X_sheep) < rang) - 1);

    % Change in population
    % Lifespan sheep
    lifespan_sheep = lifespan_sheep - 1;
    X_sheep(lifespan_sheep == 0,:) = [];
    V_sheep(lifespan_sheep == 0,:) = [];
    hunger_sheep(lifespan_sheep == 0) = [];
    lifespan_sheep(lifespan_sheep == 0) = [];

    n_sheep = length(lifespan_sheep);
    
    % Lifespan dog
    lifespan_dog = lifespan_dog - 1;
    X_dog(lifespan_dog == 0,:) = [];
    V_dog(lifespan_dog == 0,:) = [];
    hunger_dog(lifespan_dog == 0) = [];
    lifespan_dog(lifespan_dog == 0) = [];

    n_dog = length(lifespan_dog);
    
    % Starvation sheep dog
    if mod(step,5) == 0 
        hunger_sheep = hunger_sheep - 1;
        hunger_dog = hunger_dog - 1;
        
        % Sheep starve to death if hunger is -10
        X_sheep(hunger_sheep <= -10,:) = [];
        V_sheep(hunger_sheep <= -10,:) = [];
        lifespan_sheep(hunger_sheep <= -10) = [];
        hunger_sheep(hunger_sheep <= -10) = [];
        n_sheep = length(lifespan_sheep);
        
        % Dog starve to death if hunger is -10
        X_dog(hunger_dog <= -10,:) = [];
        V_dog(hunger_dog <= -10,:) = [];
        lifespan_dog(hunger_dog <= -10) = [];
        hunger_dog(hunger_dog <= -10) = [];
        n_dog = length(lifespan_dog);
    end

    % Reproduction sheep
    for i = 1:n_sheep    
        if mod(step,2) == 0 % Reproduction cycle
            for ii = 1:n_sheep
                distance_tmp = norm(X_sheep(i,:)-X_sheep(ii,:));
                if i ~= ii && distance_tmp < rang*0.1 && rand < 0.95^n_sheep ... 
                        && lifespan_sheep(i) < lifespanmax*0.8 && lifespan_sheep(ii) < lifespanmax*0.8 && hunger_sheep(i) > 0 && hunger_sheep(ii) > 0 % Do not reproduce until mature
                    hunger_sheep(i) = hunger_sheep(i) - 2;
                    hunger_sheep(ii) = hunger_sheep(ii) - 2;
                    n_sheep = n_sheep + 1;
                    X_sheep = [X_sheep; (X_sheep(i,:) + X_sheep(ii,:))/2];
                    V_sheep = [V_sheep; (V_sheep(i,:) + V_sheep(ii,:))/2];
                    lifespan_sheep = [lifespan_sheep; lifespanmax];
                    hunger_sheep = [hunger_sheep; 0];
                end
            end
        end
    end
    
    % Reproduction dog
    for j = 1:n_dog
        if mod(step,2) == 0 % Reproduction cycle
            for jj = 1:n_dog
                distance_tmp = norm(X_dog(j,:)-X_dog(jj,:));
                if j ~= jj && distance_tmp < rang*0.1 && rand < 0.99^(n_dog^1.8) ... 
                        && lifespan_dog(j) < lifespanmax*0.8 && lifespan_dog(jj) < lifespanmax*0.8  && hunger_dog(j) > 0 && hunger_dog(jj) > 0 % Do not reproduce until mature
                    hunger_dog(j) = hunger_dog(j) - 2;
                    hunger_dog(jj) = hunger_dog(jj) - 2;
                    n_dog = n_dog + 1;
                    X_dog = [X_dog; (X_dog(j,:) + X_dog(jj,:))/2];
                    V_dog = [V_dog; (V_dog(j,:) + V_dog(jj,:))/2];
                    lifespan_dog = [lifespan_dog; lifespanmax];
                    hunger_dog = [hunger_dog; 0];
                end
            end
        end
    end

    
    % Predation
    distance_sd = zeros(n_sheep,n_dog);
    distance_sg = zeros(n_sheep,n_grass);
    for i = 1:n_sheep
        for j = 1:n_dog
            distance_sd(i,j) = norm(X_sheep(i,:)-X_dog(j,:));
        end
        for k = 1:n_grass
            distance_sg(i,k) = norm(X_sheep(i,:)-X_grass(k,:));
        end
    end   

    % Predation (sheep -> grass)
    for i = 1:n_sheep
        if hunger_sheep(i) < hungermax
            X_grass(distance_sg(i,:) < rang*0.05,:) = [];
            hunger_sheep(i) = hunger_sheep(i) + length(distance_sg(:,distance_sg(i,:) < rang*0.05));
            distance_sg(:,distance_sg(i,:) < rang*0.05) = [];
            n_grass = length(X_grass(:,1));
        end
    end

    % Predation (dog -> sheep)
    for j = 1:n_dog
        if hunger_dog(j) < hungermax
            X_sheep(distance_sd(:,j) < rang*0.05,:) = [];
            V_sheep(distance_sd(:,j) < rang*0.05,:) = [];
            lifespan_sheep(distance_sd(:,j) < rang*0.05) = [];
            hunger_sheep(distance_sd(:,j) < rang*0.05) = [];
            hunger_dog(j) = hunger_dog(j) + length(distance_sd(distance_sd(:,j) < rang*0.05,:));
            distance_sd(distance_sd(:,j) < rang*0.05,:) = [];
            n_sheep = length(lifespan_sheep);
        end
    end
    
    % Update position
    X_dog = X_dog + V_dog*h;
    X_sheep = X_sheep + V_sheep*h;

    % Periodic boundary conditions
    X_dog(X_dog > rang) = X_dog(X_dog > rang) - 2*rang;
    X_dog(X_dog < -rang) = X_dog(X_dog < -rang) + 2*rang;    

    X_sheep(X_sheep > rang) = X_sheep(X_sheep > rang) - 2*rang;
    X_sheep(X_sheep < -rang) = X_sheep(X_sheep < -rang) + 2*rang;
    
    % Supply food for sheep
    if n_grass < 1000
        n_grass = n_grass + 200;
        X_grass = [X_grass; rang*(2*rand(200,2) - 1)];
    end
    
    % Save
    n_sheep_str(step) = n_sheep;
    n_dog_str(step) = n_dog;
    X_grass_str{step} = X_grass;
    X_sheep_str{step} = X_sheep;
    X_dog_str{step} = X_dog;
    
    % Termination condition
    if n_sheep == 0 || n_sheep >= 200 || n_dog == 0 || n_dog >= 30
        break
    end
end

% Result report
if step == length(t)
    disp("Sustained for " + step + " steps.")
    disp("The ecosystem was sustained!")
else
    disp("Sustained for " + step + " steps.")
end

% Plot
figure('Position',[200 200 800 400])
video = VideoWriter('ecosystem.avi','Uncompressed AVI');
open(video)
% Ecosystem simulation (left plot)
sub1 = subplot(1,2,1);
sub1.Position = [0.01 0.01 0.49 1];

s_grass = scatter(X_grass_str{1}(:,1), X_grass_str{1}(:,2), 5, [0.4660 0.6740 0.1880], 'filled', 'DisplayName', 'Grass');
hold on
s_sheep = scatter(X_sheep_str{1}(:,1), X_sheep_str{1}(:,2),20,'b','filled','DisplayName','Sheep');
s_dog = scatter(X_dog_str{1}(:,1), X_dog_str{1}(:,2),60,'r','filled','DisplayName','Dog');

axis([-rang rang -rang rang])
xticks([])
yticks([])
box on
daspect([1 1 1])
legend
hold off

% Change in population (right plot)
sub2 = subplot(1,2,2);
sub2.Position = [0.57 0.11 0.35 0.8];

yyaxis left
p_sheep = plot(0,n_sheep_str(1),'LineWidth',1,'Color','b');
xlabel('step')
ylabel('Number of non-predators (Sheep)')
axis([1 1.1 0 200])

yyaxis right
p_dog = plot(0,n_dog_str(1),'LineWidth',1,'Color','r');
xlabel('step')
ylabel('Number of predators (Dog)')
axis([1 1.1 0 30])

% Update plot
for step_plot = 2:step
    % Update scatter plot (left plot)
    s_grass.XData = X_grass_str{step_plot}(:,1);
    s_grass.YData = X_grass_str{step_plot}(:,2);

    s_sheep.XData = X_sheep_str{step_plot}(:,1);
    s_sheep.YData = X_sheep_str{step_plot}(:,2);

    s_dog.XData = X_dog_str{step_plot}(:,1);
    s_dog.YData = X_dog_str{step_plot}(:,2);
    
    % Update population plot (right plot)
    p_sheep.XData = 1:step_plot;
    p_sheep.YData = n_sheep_str(1:step_plot);
    axis([1 step_plot 0 200])

    p_dog.XData = 1:step_plot;
    p_dog.YData = n_dog_str(1:step_plot);
    axis([1 step_plot 0 30])

    frame = getframe(gcf);
    writeVideo(video,frame);
end
close(video);