/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#include <tracking/trackFindingCDC/utilities/CompositeProcessingSignalListener.h>

#include <tracking/trackFindingCDC/utilities/ReversedRange.h>
#include <iostream>

using namespace Belle2;
using namespace TrackFindingCDC;

void CompositeProcessingSignalListener::initialize()
{
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
  Super::initialize();
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
  //for (ProcessingSignalListener* psl : m_subordinaryProcessingSignalListeners) {
  //  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  //  psl->initialize();
  //}
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
}

void CompositeProcessingSignalListener::beginRun()
{
  Super::beginRun();
  for (ProcessingSignalListener* psl : m_subordinaryProcessingSignalListeners) {
    psl->beginRun();
  }
}

void CompositeProcessingSignalListener::beginEvent()
{
  Super::beginEvent();
  for (ProcessingSignalListener* psl : m_subordinaryProcessingSignalListeners) {
    psl->beginEvent();
  }
}

void CompositeProcessingSignalListener::endRun()
{
  for (ProcessingSignalListener* psl : reversedRange(m_subordinaryProcessingSignalListeners)) {
    psl->endRun();
  }
  Super::endRun();
}

void CompositeProcessingSignalListener::terminate()
{
  for (ProcessingSignalListener* psl : reversedRange(m_subordinaryProcessingSignalListeners)) {
    psl->terminate();
  }
  Super::terminate();
}

void CompositeProcessingSignalListener::addProcessingSignalListener(ProcessingSignalListener* psl)
{
  m_subordinaryProcessingSignalListeners.push_back(psl);
}

int CompositeProcessingSignalListener::getNProcessingSignalListener()
{
  return m_subordinaryProcessingSignalListeners.size();
}
