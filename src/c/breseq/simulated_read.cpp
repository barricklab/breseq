#include "libbreseq/simulated_read.h"


namespace breseq {

sim_fastq_data_t* cSimulatedFastqFactory::createFromSequence(const string &ref_sequence, const uint32_t &average_coverage, const uint32_t &read_size)
{
  cSimulatedFastqData *sim_fastq_data = new cSimulatedFastqData;
  sim_fastq_data->m_read_size = read_size;

  const size_t ref_sequence_length = ref_sequence.length();
  const uint32_t num_samples = ceil((ref_sequence_length) / read_size) * average_coverage;
  sim_fastq_data->m_sequence.reserve(num_samples * read_size);
  sim_fastq_data->m_num_N_bases = num_samples * read_size;

  // Initialize randomness
  const uint32_t &seed_value = time(NULL);
  srand(seed_value);

  sim_fastq_data->m_name = "SIMULATED-READ";

  sprintf(sim_fastq_data->m_name_plus, "SIMULATED-READ SEED:%u", seed_value);

  //As long as i is greater than 0, continue.
  //After running the loop '--' will subtract 1.
  size_t i = num_samples;
  while (i --> 0) {
    //Grab a random position based on the sequence length.
    //Potential reads outside the sequence length are forced inside.
    const size_t &start_pos = rand() % (ref_sequence_length - read_size );
    const string &segment = ref_sequence.substr(start_pos, read_size);
    sim_fastq_data->m_sequence.append(segment);
  }
  this->simulateQualityScores(sim_fastq_data);

  return sim_fastq_data;
}

void cSimulatedFastqFactory::simulateQualityScores(sim_fastq_data_t *sim_fastq_data)
{
   const uint32_t QUALITY_SCORE_OFFSET  = 33;
   const uint32_t MIN_QUALITY_SCORE  		= 0  + QUALITY_SCORE_OFFSET;
   const uint32_t MEAN_QUALITY_SCORE 		= 30 + QUALITY_SCORE_OFFSET;
   const uint32_t MAX_QUALITY_SCORE  		= 40 + QUALITY_SCORE_OFFSET;

   const size_t &length = sim_fastq_data->m_sequence.length();
   sim_fastq_data->m_qualities.resize(length);

   //! Step: Assign quality scores to each base
   cPoissonDistribution pd(MEAN_QUALITY_SCORE);
   int i = length;
   while (i-->0) {
    const char &score = char(pd.getSample(MIN_QUALITY_SCORE, MAX_QUALITY_SCORE));
    sim_fastq_data->m_qualities[i] = score;
   }
   cout << endl;
}

void cSimulatedFastqFile::write(const sim_fastq_data_t  &sim_fastq_data)
{
  const uint32_t& read_size = sim_fastq_data.m_read_size;

  ASSERT(sim_fastq_data.m_sequence.size() == sim_fastq_data.m_qualities.size(), "Error!");
  ASSERT(sim_fastq_data.m_num_N_bases % read_size == 0, "Error!");
  ASSERT(this->is_open(), "Error!");

  istringstream ss_sequence(sim_fastq_data.m_sequence);
  istringstream ss_qualities(sim_fastq_data.m_qualities);

  uint32_t current_line = 0;
  char *buffer = new char[read_size];

	const size_t &num_reads =  sim_fastq_data.m_sequence.size() / read_size;
	size_t i = num_reads;
	while (i-->0) {
		fprintf(*this, "@%s-%u\n", sim_fastq_data.m_name.c_str(), ++current_line);

    ss_sequence.get(buffer, read_size);
    fprintf(*this, "%s\n", buffer);

    fprintf(*this, "+%s\n", sim_fastq_data.m_name_plus.c_str());

    ss_qualities.get(buffer, read_size);
    fprintf(*this, "%s\n", buffer);
  }

  delete[] buffer;
}


} //end namespace breseq
