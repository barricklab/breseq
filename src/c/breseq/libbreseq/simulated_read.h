#ifndef SIMULATED_READ_H
#define SIMULATED_READ_H

#include "fastq.h"
#include "common.h"

namespace breseq {

class cSimulatedFastqFactory
{
public:
  struct cSimulatedFastqData : public cFastqSequence{};

  cSimulatedFastqFactory(){}

  cSimulatedFastqData* createFromSequence(const string &sequence, const uint32_t &average_coverage);

private:
  void simulateQualityScores(cSimulatedFastqData *sim_fastq_data);
};

typedef cSimulatedFastqFactory sim_fastq_factory_t;
typedef cSimulatedFastqFactory::cSimulatedFastqData sim_fastq_data_t;

class cSimulatedFastqFile : public ofstream
{
  public:
    cSimulatedFastqFile(const string &file_name)
      : ofstream(file_name.c_str()) {}

    void write(const sim_fastq_data_t &sim_fastq_data);
};

} //end namespace breseq

#endif // SIMULATED_READ_H
