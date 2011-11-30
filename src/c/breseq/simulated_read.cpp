#include "libbreseq/simulated_read.h"


namespace breseq {

sim_fastq_data_t* cSimulatedFastqFactory::createFromSequence(const string &sequence, const uint32_t &average_coverage)
{
  cSimulatedFastqData *sim_fastq_data = new cSimulatedFastqData;

  const size_t &sequence_length = sequence.length();
  uint32_t n_samples = floor((sequence_length - 1) / 36) * average_coverage;
  sim_fastq_data->m_sequence.resize(n_samples * 36);
  sim_fastq_data->m_num_N_bases = n_samples * 36;

  //Lazy fix for handling ends for now
  string base_sequence = sequence + string(35, 'N');

  srand(time(NULL));

  string::const_iterator it_base;
  string::iterator it_sim = sim_fastq_data->m_sequence.begin();

  uint32_t sample_pos;

  while (n_samples --> 0) {
    it_base = base_sequence.begin();

    sample_pos = rand() % sequence_length;
    advance(it_base, sample_pos);

    copy(it_base, it_base + 35, it_sim);
    advance(it_sim, 35);
  }
  this->simulateQualityScores(sim_fastq_data);

  return sim_fastq_data;
}

void cSimulatedFastqFactory::simulateQualityScores(cSimulatedFastqData *sim_fastq_data)
{
  const size_t &length = sim_fastq_data->m_sequence.length();
  sim_fastq_data->m_qualities.resize(length);
  sim_fastq_data->m_qualities.assign(length, '#');
}


void cSimulatedFastqFile::write(const cSimulatedFastqFactory::cSimulatedFastqData &sim_fastq_data)
{
  ASSERT(sim_fastq_data.m_sequence.size() == sim_fastq_data.m_qualities.size(), "Error!");
  ASSERT(sim_fastq_data.m_num_N_bases % 36 == 0, "Error!");
  ASSERT(this->is_open(), "Error!");

  istringstream ss_sequence(sim_fastq_data.m_sequence);
  istringstream ss_qualities(sim_fastq_data.m_qualities);

  uint32_t current_line = 0;

  char *buffer = new char[36];
  buffer[35] = '\n';

  const uint32_t &num_reads = sim_fastq_data.m_num_N_bases / 36;
  for (int i = 0; i < num_reads; i++) {
    (*this) <<  "@READ-" << ++current_line << endl;

    ss_sequence.read(buffer, 35);
    ofstream::write(buffer, 36);

    ofstream::write("+\n", 2);

    ss_qualities.read(buffer, 35);
    ofstream::write(buffer, 36);
  }

  delete[] buffer;


}


} //end namespace breseq
