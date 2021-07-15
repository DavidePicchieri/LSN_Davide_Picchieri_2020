#include "genetic.h"

using namespace std;



Specimen :: Specimen(int ncities, Random * rnd){
  for(int i=0; i<ncities; i++) m_cities.push_back(i);
  m_rnd=rnd;
}

Specimen :: ~Specimen(){}


Specimen& Specimen :: operator=( Specimen& a) // Overload
{
  for(int i=0; i<m_cities.size(); i++){
    m_cities[i]=a.GetCities()[i];
  }
  m_distance = a.m_distance;
  m_rnd = a.m_rnd;
  m_gen = a.m_gen;
  
  return *this;
}

void Specimen :: Check (){ //Check on the trip, if there are repetitions
  vector<int> sorted;
  for(int i=0; i<m_cities.size(); i++) sorted.push_back(m_cities[i]);
  sort(sorted.begin(), sorted.end());
  for(int i=0; i< m_cities.size()-1; i++){
    if(sorted[i]==sorted[i+1]){
      cout << "Two cities are the same!" << endl;
      abort();
    }
  }
}


void Specimen :: Print(string filename){ //Print single trip
  ofstream Output;
  Output.open(filename, ios::app);
  for(int i=0; i<m_cities.size(); i++){
    Output << m_cities[i] << endl;
  }
  Output << endl << endl;
}

int Specimen :: Pbc(int i){ //Periodic boundary conditions without the 0th element
  if(i >  m_cities.size()){
    i = i - m_cities.size()+1; 
  }
  else if(i == m_cities.size()){
    i = i - m_cities.size()+1;
  } 
  else if(i == 0){
    i = i + m_cities.size()-1;
  } 
  else if(i < 0){  
    i = i + m_cities.size()-1;
  }
  return i;
}

void Specimen :: Swap(int i){ //Swap two cities
  if(i>=m_cities.size()) {
    cout << "Error i >= N" << endl;
    abort();
  }
  
  if(i==0) {
    cout << "Error i=0" << endl;
    abort();
  }
  int appo=m_cities[i];
  m_cities[i]=m_cities[Pbc(i+1)];
  m_cities[Pbc(i+1)]=appo;
}

void Specimen :: Shift(int i, int groupsize, int shift){ //Shift a group of cities
  if(groupsize>=m_cities.size()-1) {
    cout << "Error groupsize>=N-1" << endl;
    abort();
  }
  if(i>=m_cities.size() or i==0) {
    cout << "Error i>=N or i=0" << endl;
    abort();
  }
  if(shift>=m_cities.size()-1) {
    cout << "Error shift>=N-1" << endl;
    abort();
  }
  
  if(shift+groupsize>=m_cities.size()-1) {
    cout << "Error shift+groupsize>=N" << endl;
    abort();
  }
  
  int appo1[shift+groupsize];
  int appo2[shift+groupsize];
  
  for(int j=0; j<groupsize; j++) appo1[j]=m_cities[Pbc(i+j)];
  for(int j=0; j<shift; j++) appo2[j]=m_cities[Pbc(Pbc(i+groupsize)+j)];
  for(int j=0; j<shift; j++) m_cities[Pbc(i+j)]=appo2[j];
  for(int j=0; j<groupsize; j++) m_cities[Pbc(Pbc(i+shift)+j)]=appo1[j];
  
}

void Specimen :: SwapGroup(int i, int groupsize){ //Swap a group of cities 
  if(groupsize>=m_cities.size()/2) {
    cout << "Error groupsize>=N/2" << endl;
    abort();
  }
  if(i>=m_cities.size() or i==0) {
    cout << "Error i>=N or i=0" << endl;
    abort();
  }
  
  int appo1[groupsize];
  int appo2[groupsize];
  
  for(int j=0; j<groupsize; j++){
    appo1[j]=m_cities[Pbc(i+j)];
    appo2[j]=m_cities[Pbc(Pbc(i+groupsize)+j)];
    m_cities[Pbc(Pbc(i+groupsize)+j)]=appo1[j];
    m_cities[Pbc(i+j)]=appo2[j];
  }
}

void Specimen :: Inversion(int i, int groupsize){ //Invert a group of cities
  if(groupsize>=m_cities.size()-1) {
    cout << "Error groupsize>=N-1" << endl;
    abort();
  }
  if(i>=m_cities.size() or i==0) {
    cout << "Error i>=N or i=0" << endl;
    abort();
  }
  
  vector<int> appo;
  for(int j=0; j<groupsize; j++) appo.push_back(m_cities[Pbc(i+j)]);
  std::reverse(appo.begin(), appo.end());
  for(int j=0; j<groupsize; j++) m_cities[Pbc(i+j)]=appo[j];
}

///////////////////SPECIES////////////////////////////////

Species :: Species(int popsize, Random * rnd, int ncities){
    m_popsize = popsize;
    m_rnd = rnd;
    m_ncities = ncities;
    for(int i=0; i<m_popsize; i++){
        Specimen appo(m_ncities, m_rnd);
        m_pop.push_back(appo);
        m_popnew.push_back(appo);
    }
    m_ngen = 0;
}
Species :: ~Species(){}

void Species :: Print(int index, string filename){ //Print a specific trip
  m_pop[index].Print(filename);
}
void Species :: PrintBest(string filename){ //Print the coords of the best trip
  ofstream Output;
  Output.open(filename);  
  for(int i=0; i<m_ncities; i++){
    Output << m_besttrip[i] << "\t" << m_pos[0][m_besttrip[i]] << "\t" << m_pos[1][m_besttrip[i]] << endl;
  }
  Output << m_besttrip.size() << "\t" << m_pos[0][0] << "\t" << m_pos[1][0] << endl;
  Output.close();
}

void Species :: PrintDistance(string filename){ //Overloaded to print whole population or single individual
  ofstream Output;
  Output.open(filename);
  if (!Output.is_open()){
      abort;
  }
  for (int i=0; i<m_popsize; i++){
      Output << m_pop[i].GetDistance() << endl;
  }
  Output.close();
}
void Species :: PrintDistance(int index, string filename){
  ofstream Output;
  Output.open(filename, ios::app);
  Output << m_pop[index].GetDistance() << endl;
  Output.close();
}
void Species :: PrintAve(int nspecimen, string filename){//Print the average distance of the best nspecimen trips of the population
  double ave = 0;
  if (nspecimen > m_popsize){
    cout << "Insert < " << m_popsize << " as argument of the function PrintAve" << endl;
    abort;
  }
  for (int i=0; i<nspecimen; i++){
    ave += m_pop[i].GetDistance();
  }
  ave /= double(nspecimen);
  ofstream Output;
  Output.open(filename, ios::app);
  Output << m_ngen << "\t" << ave << endl;
  Output.close();
}
void Species :: PrintBestSA(string filename){ 
  ofstream Output;
  Output.open(filename);
  
  for(int i=0; i<m_ncities; i++){
    Output << m_pop[0].GetCities()[i] << setw(12) << m_pos[0][m_pop[0].GetCities()[i]] << setw(12) << m_pos[1][m_pop[0].GetCities()[i]] << endl;
  }
  Output << m_pop[0].GetCities().size() << setw(12) << m_pos[0][0] << setw(12) << m_pos[1][0] << endl;
}

void Species::Copy(int i, int j){ //Copy two individuals of the population
    m_pop[i]=m_pop[j];
}  
void Species :: SetPosCircle(){ //Initialize on a circle
  vector<double> x, y;
  double angle, r=1;
  x.push_back(0); //initial position always set in (0,1)
  y.push_back(1);

  for (int i=1; i<m_ncities; i++){
      angle = m_rnd->Rannyu(0,2*M_PI);
      x.push_back(r*cos(angle));
      y.push_back(r*sin(angle));
  }
    
  m_pos.push_back(x);
  m_pos.push_back(y);

  for(int i=0; i<1000; i++){ //Swapping many times the cities in the different specimen to create the initial population
    for(int j=0; j<m_popsize; j++){
      m_pop[j].Swap(int(m_rnd->Rannyu(1, m_ncities)));
    }
  }
}

void Species :: SetPosSquare(){ //Initialize on a square
  double side = 1;
  vector<double> x, y;
  for (int i=0; i<m_ncities; i++){
      x.push_back(m_rnd->Rannyu(-side, side));
      y.push_back(m_rnd->Rannyu(-side, side));
  }

  m_pos.push_back(x);
  m_pos.push_back(y);

  for(int i=0; i<1000; i++){ //Swapping many times the cities in the different specimen to create the initial population
    for(int j=0; j<m_popsize; j++){
      m_pop[j].Swap(int(m_rnd->Rannyu(1, m_ncities)));
    }
  }
}

void Species :: Measure(){ //Measure the distances of the trips in the population
  int appo1, appo2;
  double distance, dx, dy;

  for(int i=0; i<m_popsize; i++){
    distance=0;

    for(int j=0; j<m_ncities-1; j++){
      appo1 = m_pop[i].GetCities()[j];
      appo2 = m_pop[i].GetCities()[j+1];
      dx = m_pos[0][appo1] - m_pos[0][appo2];
      dy = m_pos[1][appo1] - m_pos[1][appo2];
      distance += sqrt(dx*dx+dy*dy);
    }

    appo1 = m_pop[i].GetCities()[m_ncities-1];
    appo2 = m_pop[i].GetCities()[0];
    dx = m_pos[0][appo1] - m_pos[0][appo2];
    dy = m_pos[1][appo1] - m_pos[1][appo2];
    distance += sqrt(dx*dx+dy*dy);
    
    m_pop[i].SetDistance(distance);
  }

}
void Species :: Sort(){ //Sort the trips so that the shortest is in position [0]
  int index;
  m_besttrip = m_pop[0].GetCities();
  for (int i = 0; i < m_popsize-1; i++){
    index = i;
    for (int j = i+1; j < m_popsize; j++){
      if (m_pop[j].GetDistance() < m_pop[index].GetDistance()) 
        index = j;
    }
    Specimen appo = m_pop[i];
    m_pop[i] = m_pop[index];
    m_pop[index] = appo;
  }
  if(m_pop[0].GetDistance()<m_mindist){
    m_besttrip = m_pop[0].GetCities();
    m_mindist = m_pop[0].GetDistance();
    m_bestgen = m_ngen;
  }
}
void Species :: Selection(double powerlaw){ //Select two trips to go to the next generation by using a power law and their position in the vector which is related to their length
  int index1, index2;
  for(int i=0; i<m_popsize; i=i+2){
    
    index1=int(m_popsize*pow(m_rnd->Rannyu(),powerlaw));
    do {
      index2=int(m_popsize*pow(m_rnd->Rannyu(),powerlaw));
    }
    while(index2==index1);
    
    m_popnew[i]=m_pop[index1];
    m_popnew[i+1]=m_pop[index2];
  }
}
void Species :: Crossover(int cut, double chance){ //Puts the trips in the new generation back in the main one with a chance of a crossover happening
  int count1, count2;
  for(int i=0; i<m_popsize; i=i+2){
      m_pop[i]=m_popnew[i];
      m_pop[i+1]=m_popnew[i+1];
      if(m_rnd->Rannyu() < chance){
        
        int missing1[m_ncities-cut];
        int missing2[m_ncities-cut];
        for(int j=0; j<m_ncities-cut; j++){
          missing1[j]=m_popnew[i].GetCities()[j+cut];
          missing2[j]=m_popnew[i+1].GetCities()[j+cut];
        }    
        count1=cut;
        count2=cut;
        for(int j=1; j<m_ncities; j++){ 
          for(int k=0; k<m_ncities-cut; k++){
            if(m_popnew[i+1].GetCities()[j] == missing1[k]){
              m_pop[i].SetCity(count1, missing1[k]);
              count1++;
            }
            if(m_popnew[i].GetCities()[j] == missing2[k]){
              m_pop[i+1].SetCity(count2, missing2[k]);
              count2++;
            }
          }
        }
        for(int k=0; k<m_popsize; k++) m_pop[k].Check();
      }
    }
  m_ngen=m_ngen+1;

}
void Species :: Mutation(double chance){ //Calls the mutation methods for the trips with a certain chance
  int groupsize, shift;

  for(int i=0; i<m_popsize; i++){  //Swap
    if(m_rnd->Rannyu()<chance){
      m_pop[i].Swap(int(m_rnd->Rannyu(1, m_ncities)));
    }
  }

  for(int i=0; i<m_popsize; i++){  //Shift
      if(m_rnd->Rannyu()<chance){
        do{
          groupsize = int(m_rnd->Rannyu(1,m_ncities-1));
          shift = int(m_rnd->Rannyu(1,m_ncities-1));
        }while(groupsize + shift >=(m_ncities-1));

        m_pop[i].Shift(int(m_rnd->Rannyu(1, m_ncities)), groupsize, shift);
    }
   }

  for(int i=0; i<m_popsize; i++){  //SwapGroup
      if(m_rnd->Rannyu()<chance){
      groupsize = int(m_rnd->Rannyu(1, m_ncities/2));
      m_pop[i].SwapGroup(int(m_rnd->Rannyu(1, m_ncities)), groupsize);
    }
  }
  
  for(int i=0; i<m_popsize; i++){  //Inversion
      if(m_rnd->Rannyu()<chance){
      groupsize = int(m_rnd->Rannyu(1, m_ncities-1));
      m_pop[i].Inversion(int(m_rnd->Rannyu(1, m_ncities)), groupsize);
    }
  }

  for(int i=0; i<m_popsize; i++) m_pop[i].Check();
}
  ////////methods used in exercise 10.1////

void Species :: SetCoeff(double coeff){ 
    m_ann_coeff = coeff;
}
void Species :: SetTemp(double temp){
  m_temp = temp;
}
void Species :: Annealing(){
  m_temp *= m_ann_coeff;
}
void Species :: Accept(){
  if(m_rnd->Rannyu()<this->GetAlpha()){
    m_pop[0]=m_pop[1];
  }else{
    m_pop[1]=m_pop[0];
  }
  m_ngen=m_ngen+1;
}
double Species :: GetAlpha(){
  return fmin(1., exp(-1./m_temp*(m_pop[1].GetDistance()-m_pop[0].GetDistance())));
}
void Species :: Rearrange(){
  int j = m_rnd->Rannyu(1, m_ncities);
  m_pop[1].Swap(j);

  j=int(m_rnd->Rannyu(1, m_ncities));
  int m = m_rnd->Rannyu(1, m_ncities-1);
  m_pop[1].Inversion(j, m);
  
  m_pop[1].Check();
}

