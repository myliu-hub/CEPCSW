/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#include "CompositeProcessingSignalListener.h"

#include "ReversedRange.h"
#include <iostream>

using namespace Belle2;
using namespace TrackFindingCDC;

void CompositeProcessingSignalListener::initialize()
{
  Super::initialize();
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
