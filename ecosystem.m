clear
close all

%各エージェントの初期個体数
n_grass = 200;
n_sheep = 100;
n_dog = 10;

%ステップ幅，時間
h = 0.01;
t = 0:h:30;

%初期状態（位置・速度）
rang = 10; 

X_grass = rang*(2*rand(n_grass,2) - 1); %grass 初期位置
X_sheep = 2*(2*rand(n_sheep,2) - 1); %sheep 初期位置
X_dog = rang*rand(n_dog,2) - rang; %dog 初期位置

V_sheep = 2*rand(n_sheep,2) - 1; %sheep 初期速度
V_dog = 2*rand(n_dog,2) - 1; %dog 初期速度

Vs = 3*(2*rand(1,2) - 1); %sheepの進行速度

%寿命設定
lifespanmax = 300; %寿命最大
lifespan_sheep = lifespanmax*ones(n_sheep,1);
lifespan_dog = lifespanmax*ones(n_dog,1);

%空腹度設定
hungermax = 20; %満腹状態
hunger_sheep = hungermax*ones(n_sheep,1);
hunger_dog = randi([1 10],n_dog,1);

%各ゲイン
k_sheep_av = 0;
k_sheep_rp = 1.2;
k_sheep_at = 0.2;

k_dog_at = 0.1;

%各値の保存用
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

%シミュレーション
disp("Give me a second...")
for step = 2:length(t)
    %速度更新
    V_sheep = V_sheep - k_sheep_av*(V_sheep - sum(V_sheep)/n_sheep); %Void(整列)
    V_sheep = V_sheep - k_sheep_av*(V_sheep - sum(X_sheep)/n_sheep); %Void(結合)

    if mod(step,30) == 0 %進行方向(サイクル30)
        Vs = 2*rand(1,2) - 1;
    end

    V_sheep = V_sheep + 3*(2*rand(size(V_sheep)) - 1) + ones(size(V_sheep)).*Vs; %ランダム＋進行方向
    V_dog = V_dog + 10*(2*rand(size(V_dog)) - 1); %ランダム

    for i = 1:n_sheep
        for j = 1:n_dog
            distance_tmp = norm(X_sheep(i,:)-X_dog(j,:));
            if distance_tmp < rang && hunger_dog(j) <= 0 %飢餓状態
                V_dog(j,:) = V_dog(j,:) - k_dog_at*(X_dog(j,:) - X_sheep(i,:))/distance_tmp;
            end
            if distance_tmp < rang*0.6 %回避
                V_sheep(i,:) = V_sheep(i,:) + k_sheep_rp*(X_sheep(i,:) - X_dog(j,:))/distance_tmp;
            end
        end
        for ii = 1:n_sheep
            distance_tmp = norm(X_sheep(i,:)-X_sheep(ii,:));
            if i ~= ii && distance_tmp > rang*0.5 && distance_tmp < rang*1.2 %同種に近づく
                V_sheep(i,:) = V_sheep(i,:) - k_sheep_at*(X_sheep(i,:) - X_sheep(ii,:))/distance_tmp;
            end
        end
        
    end
    V_sheep = V_sheep*0.9 ; %減速
    V_dog = V_dog*0.9;

    %壁境界条件
%     V_dog = V_dog.*(2*(abs(X_dog) < rang) - 1);
%     V_sheep = V_sheep.*(2*(abs(X_sheep) < rang) - 1);

    %各個体数の変動
    %寿命sheep
    lifespan_sheep = lifespan_sheep - 1;
    X_sheep(lifespan_sheep == 0,:) = [];
    V_sheep(lifespan_sheep == 0,:) = [];
    hunger_sheep(lifespan_sheep == 0) = [];
    lifespan_sheep(lifespan_sheep == 0) = [];

    n_sheep = length(lifespan_sheep);
    
    %寿命dog
    lifespan_dog = lifespan_dog - 1;
    X_dog(lifespan_dog == 0,:) = [];
    V_dog(lifespan_dog == 0,:) = [];
    hunger_dog(lifespan_dog == 0) = [];
    lifespan_dog(lifespan_dog == 0) = [];

    n_dog = length(lifespan_dog);
    
    %餓死sheep dog
    if mod(step,5) == 0 
        hunger_sheep = hunger_sheep - 1;
        hunger_dog = hunger_dog - 1;
        
        %sheep 空腹度が－10で餓死
        X_sheep(hunger_sheep <= -10,:) = [];
        V_sheep(hunger_sheep <= -10,:) = [];
        lifespan_sheep(hunger_sheep <= -10) = [];
        hunger_sheep(hunger_sheep <= -10) = [];
        n_sheep = length(lifespan_sheep);
        
        %dog 空腹度が－10で餓死
        X_dog(hunger_dog <= -10,:) = [];
        V_dog(hunger_dog <= -10,:) = [];
        lifespan_dog(hunger_dog <= -10) = [];
        hunger_dog(hunger_dog <= -10) = [];
        n_dog = length(lifespan_dog);
    end

    %繁殖sheep
    for i = 1:n_sheep    
        if mod(step,5) == 0 %繁殖周期
            for ii = 1:n_sheep
                distance_tmp = norm(X_sheep(i,:)-X_sheep(ii,:));
                if i ~= ii && distance_tmp < rang*0.1 && rand < 0.95^n_sheep ... 
                        && lifespan_sheep(i) < lifespanmax*0.8 && lifespan_sheep(ii) < lifespanmax*0.8 %成熟するまで繁殖しない
                    n_sheep = n_sheep + 1;
                    X_sheep = [X_sheep; (X_sheep(i,:) + X_sheep(ii,:))/2];
                    V_sheep = [V_sheep; (V_sheep(i,:) + V_sheep(ii,:))/2];
                    lifespan_sheep = [lifespan_sheep; lifespanmax];
                    hunger_sheep = [hunger_sheep; hungermax];
                end
            end
        end
    end
    
    %繁殖dog
    for j = 1:n_dog
        if mod(step,3) == 0 %繁殖周期
            for jj = 1:n_dog
                distance_tmp = norm(X_dog(j,:)-X_dog(jj,:));
                if j ~= jj && distance_tmp < rang*0.1 && rand < 0.9^n_dog ... 
                        && lifespan_dog(j) < lifespanmax*0.8 && lifespan_dog(jj) < lifespanmax*0.8 %成熟するまで繁殖しない
                    n_dog = n_dog + 1;
                    X_dog = [X_dog; (X_dog(j,:) + X_dog(jj,:))/2];
                    V_dog = [V_dog; (V_dog(j,:) + V_dog(jj,:))/2];
                    lifespan_dog = [lifespan_dog; lifespanmax];
                    hunger_dog = [hunger_dog; 0];
                end
            end
        end
    end

    
    %捕食
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

    %捕食(sheep　→　grass)
    for i = 1:n_sheep
        X_grass(distance_sg(i,:) < rang*0.05,:) = [];
        hunger_sheep(i) = hunger_sheep(i) + length(distance_sg(:,distance_sg(i,:) < rang*0.05));
        distance_sg(:,distance_sg(i,:) < rang*0.05) = [];
        n_grass = length(X_grass(:,1));
    end

    %捕食(dog → sheep)
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
    
    %位置更新
    X_dog = X_dog + V_dog*h;
    X_sheep = X_sheep + V_sheep*h;

    %周期的境界条件
    X_dog(X_dog > rang) = X_dog(X_dog > rang) - 2*rang;
    X_dog(X_dog < -rang) = X_dog(X_dog < -rang) + 2*rang;    

    X_sheep(X_sheep > rang) = X_sheep(X_sheep > rang) - 2*rang;
    X_sheep(X_sheep < -rang) = X_sheep(X_sheep < -rang) + 2*rang;
    
    %sheepの餌供給
    if n_grass < 600
        n_grass = n_grass + 100;
        X_grass = [X_grass; rang*(2*rand(100,2) - 1)];
    end
    
    %保存
    n_sheep_str(step) = n_sheep;
    n_dog_str(step) = n_dog;
    X_grass_str{step} = X_grass;
    X_sheep_str{step} = X_sheep;
    X_dog_str{step} = X_dog;
    
    %終了条件
    if n_sheep == 0 || n_sheep >= 200 || n_dog == 0 || n_dog >= 20
        break
    end
end

%結果報告
if step == length(t)
    disp("Sustained for " + step + " steps.")
    disp("The ecosystem was sustained!")
else
    disp("Sustained for " + step + " steps.")
end

%プロット
figure('Position',[200 200 800 400])
video = VideoWriter('ecosystem.mp4','MPEG-4');
open(video);

%生態系シミュレーション(左図)
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

%個体数の推移(右図)
sub2 = subplot(1,2,2);
sub2.Position = [0.57 0.11 0.35 0.8];

yyaxis left
p_sheep = plot(0,n_sheep_str(1),'LineWidth',1,'Color','b');
xlabel('step')
ylabel('Number of non-predators(Sheep)')
axis([1 1.1 0 200])

yyaxis right
p_dog = plot(0,n_dog_str(1),'LineWidth',1,'Color','r');
xlabel('step')
ylabel('Number of predators(Dog)')
axis([1 1.1 0 20])


%プロット更新
for step_plot = 2:step
    %点群の更新(左図)
    s_grass.XData = X_grass_str{step_plot}(:,1);
    s_grass.YData = X_grass_str{step_plot}(:,2);

    s_sheep.XData = X_sheep_str{step_plot}(:,1);
    s_sheep.YData = X_sheep_str{step_plot}(:,2);

    s_dog.XData = X_dog_str{step_plot}(:,1);
    s_dog.YData = X_dog_str{step_plot}(:,2);
    
    %個体数推移の更新(右図)
    p_sheep.XData = 1:step_plot;
    p_sheep.YData = n_sheep_str(1:step_plot);
    axis([1 step_plot 0 200])

    p_dog.XData = 1:step_plot;
    p_dog.YData = n_dog_str(1:step_plot);
    axis([1 step_plot 0 20])

    frame = getframe(gcf);
    writeVideo(video,frame);  
end
close(video);