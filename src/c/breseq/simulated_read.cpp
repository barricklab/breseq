#include "libbreseq/simulated_read.h"


namespace breseq {

sim_fastq_data_t* cSimulatedFastqFactory::createFromSequence(const string &sequence, const uint32_t &average_coverage)
{
  cSimulatedFastqData *sim_fastq_data = new cSimulatedFastqData;

  const size_t &sequence_length = sequence.length();
  uint32_t n_samples = ceil((sequence_length) / mReadLength) * average_coverage;
  sim_fastq_data->m_sequence.resize(n_samples * mReadLength);
  sim_fastq_data->m_num_N_bases = n_samples * mReadLength;

  // Initialize randomness
  srand(time(NULL));

  string::const_iterator it_base;
  string::iterator it_sim = sim_fastq_data->m_sequence.begin();

  uint32_t sample_pos;

  //As long as n_samples is greater than 0, continue.
  //After running the loop '--' will subtract 1.
  while (n_samples-- > 0) {
    //String iterator is set to the first position.
    //This position is position 1.
    it_base = sequence.begin();

    //Grab a random position based on the sequence length.
    //Potential reads outside the sequence length are forced inside.
    sample_pos = rand() % (sequence_length + (mReadLength - 1));
    sample_pos = max(sample_pos, mReadLength - 1);
    sample_pos = min(sample_pos, (uint32_t)(sequence_length));
    sample_pos -= (mReadLength - 1);
    
    //Advance the iterator to the position we generated.
    advance(it_base, sample_pos);

    copy(it_base, it_base + (mReadLength - 1), it_sim);
    advance(it_sim, (mReadLength - 1));
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

void cSimulatedFastqFile::write(const cSimulatedFastqFactory::cSimulatedFastqData &sim_fastq_data, uint32_t uReadLength)
{
  ASSERT(sim_fastq_data.m_sequence.size() == sim_fastq_data.m_qualities.size(), "Error!");
  ASSERT(sim_fastq_data.m_num_N_bases % uReadLength == 0, "Error!");
  ASSERT(this->is_open(), "Error!");

  istringstream ss_sequence(sim_fastq_data.m_sequence);
  istringstream ss_qualities(sim_fastq_data.m_qualities);

  uint32_t current_line = 0;

  char *buffer = new char[uReadLength];
  buffer[(uReadLength - 1)] = '\n';

  const uint32_t &num_reads = sim_fastq_data.m_num_N_bases / uReadLength;
  for (uint32_t i = 0; i < num_reads; i++) {
    (*this) <<  "@READ-" << ++current_line << endl;

    ss_sequence.read(buffer, (uReadLength - 1));
    ofstream::write(buffer, uReadLength);

    ofstream::write("+\n", 2);

    ss_qualities.read(buffer, (uReadLength - 1));
    ofstream::write(buffer, uReadLength);
  }

  delete[] buffer;
}


} //end namespace breseq
