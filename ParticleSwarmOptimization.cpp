#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<random>
#include<algorithm>


using namespace std;

struct City{
    int id;
    double x;
    double y;
};

struct Swap{
    int i;
    int j;
};

struct Particle {
    vector<int> position;
    Swap velocity;
    vector<int> pbest;
    double pbest_fitness;
    vector<int> informants;
};

double fitness(const vector<vector<double>>& distance_matrix, const vector<int>& tour){
    double total_distance = 0.0;

    for (int i = 0; i < tour.size() - 1; i++){
        total_distance += distance_matrix[tour[i]][tour[i+1]];
    }
    total_distance += distance_matrix[tour.back()][tour[0]];

    return 1 / total_distance;
}

vector<int> buildTour(const vector<City>& cities){
    vector<int> tour(cities.size());
    for (int i = 0; i < cities.size(); i++){
        tour[i] = i;
    }
    shuffle(tour.begin(), tour.end(), mt19937(random_device()()));
    return tour;
}

Swap randomSwap(int n, mt19937& gen){
    uniform_int_distribution<int> dist(0, n - 1);

    int i = dist(gen);
    int j = dist(gen);

    while (i == j) {
        j = dist(gen);
    }
    return {i, j};
}

void applySwap(vector<int>& position, const Swap& velocity){
    swap(position[velocity.i], position[velocity.j]);
}

Swap findSwap(const vector<int>& position, const vector<int>& target){
    for (int i = 0; i < position.size(); i++){
        if (position[i] != target[i]){
            for (int j = i + 1; j < position.size(); j++){
                if (position[j] == target[i]){
                    return {i, j};
                }
            }
        }
    }
    return {-1, -1};
}

vector<int> chooseInformants(int self_index, int swarm_size, int k, mt19937& gen){
    vector<int> pool;
    for (int i = 0; i < swarm_size; i++){
        if (i != self_index) pool.push_back(i);
    }
    shuffle(pool.begin(), pool.end(), gen);
    if ((int)pool.size() > k) pool.resize(k);
    return pool;
}

vector<int> bestInformantPosition(const vector<Particle>& particles, const vector<int>& informants){
    double best_fitness = -1.0;
    int best_idx = informants[0];
    for (int idx : informants){
        if (particles[idx].pbest_fitness > best_fitness){
            best_fitness = particles[idx].pbest_fitness;
            best_idx = idx;
        }
    }
    return particles[best_idx].pbest;
}

double PSO(const vector<vector<double>>& distance_matrix,
           vector<Particle>& particles,
           int max_iterations){

    double alpha = 0.5;
    double beta = 0.5;
    double gamma = 0.5;
    double delta = 0.5;
    double epsilon = 0.1;

    int informant_count = 4; 

    vector<int> gbest;
    double gbest_fitness = 0.0;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);

    for (int i = 0; i < particles.size(); i++){
        if (particles[i].pbest_fitness > gbest_fitness){
            gbest_fitness = particles[i].pbest_fitness;
            gbest = particles[i].pbest;
        }
    }
    
    for (int k = 0; k < max_iterations; k++){
        for (size_t idx = 0; idx < particles.size(); idx++){
            Particle& particle = particles[idx];

            particle.informants = chooseInformants((int)idx, (int)particles.size(), informant_count, gen);

            double b = dist(gen);
            double c = dist(gen);
            double a = dist(gen);
            double d = dist(gen);
            double e = dist(gen);

            if (e < epsilon){
                Swap random_swap = randomSwap(particle.position.size(), gen);
                applySwap(particle.position, random_swap);
            }

            if (d < delta){
                vector<int> informant_best = bestInformantPosition(particles, particle.informants);
                Swap informant_swap = findSwap(particle.position, informant_best);
                if (informant_swap.i != -1){
                    applySwap(particle.position, informant_swap);
                }
            }

            if (a < alpha){
                if (particle.velocity.i != -1){
                    applySwap(particle.position, particle.velocity);
                }
            }

            if (b < beta){
                Swap pbest_swap = findSwap(particle.position, particle.pbest);
                if (pbest_swap.i != -1){
                    applySwap(particle.position, pbest_swap);
                    particle.velocity = pbest_swap;
                }
            }
            if (c < gamma){
                Swap gbest_swap = findSwap(particle.position, gbest);
                if (gbest_swap.i != -1){
                    applySwap(particle.position, gbest_swap);
                }
            }

            double current_fitness = fitness(distance_matrix, particle.position);

            if (current_fitness > particle.pbest_fitness){
                particle.pbest_fitness = current_fitness;
                particle.pbest = particle.position;
            }
            if (particle.pbest_fitness > gbest_fitness){
                gbest_fitness = particle.pbest_fitness;
                gbest = particle.pbest;
            }

        }
        cout << "Iteration " << k + 1 << " | Best fitness: " << gbest_fitness << endl;
    }

    return gbest_fitness;
}

int main(){
    ifstream instance("data/berlin52.tsp");

    if (!instance.is_open()){
        cout << "Error opening file check extension" << endl;
        return 1;
    }

    string line;
    vector<City> cities;
    bool done = false;

    while(getline(instance, line)){
        if (line != "NODE_COORD_SECTION"){
            continue;
        }else{
            while(getline(instance, line)){
                if (line == "EOF"){
                    done = true;
                    break;
                }
                City city;
                sscanf(line.c_str(), "%d %lf %lf", &city.id, &city.x, &city.y);
                cities.push_back(city);
            }
        }
        if (done) break;
    }

    vector<vector<double>> distance_matrix(cities.size(), vector<double>(cities.size()));

    for (int i = 0; i < cities.size(); i++){
        for (int j = 0; j < cities.size(); j++){
            double dx = cities[i].x - cities[j].x;
            double dy = cities[i].y - cities[j].y;
            double d = sqrt(dx * dx + dy * dy);
            distance_matrix[i][j] = d;
            distance_matrix[j][i] = d;
        }
    }

    random_device rd;
    mt19937 gen(rd());

    vector<Particle> particles;
    int population_size = 30;

    for (int i = 0; i < population_size; i++){
        Particle particle;

        particle.position = buildTour(cities);
        particle.pbest = particle.position;
        particle.pbest_fitness = fitness(distance_matrix, particle.position);
        particle.velocity = randomSwap(cities.size(), gen);

        particles.push_back(particle);
    }
    double result = PSO(distance_matrix, particles, 100);

    cout << "Best fitness found: " << result << endl;

    return 0;
}
