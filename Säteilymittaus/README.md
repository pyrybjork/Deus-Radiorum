# Säteilymittaus

Ideana oli verrata meidän mittauksia Segeyltä saatuun malliin, joka esittää täälläpäin havaittavaa säteilyä korkeuteen verrattuna. Ongelmana on kuitenkin, että ilmeisesti meidän muovi näyttää vain myonisäteilyä ja tuo malli näyttää kokonaissäteilyn. 

Sergey-kansiossa on sen lähettämät scriptit, mutta en saanut niitä toimimaan. Kysyin Sergeyltä, unohtiko se lähettää tuon _common kirjaston vai enkä vaan osannut.

## Tietoa Sergeyn mallista

This radiation originates from so-called "showers" of different secondary particles (protons, neutrons, electrons, gamma-rays) which are the product of interaction between cosmic rays and nuclei of the atmosphere.

My colleagues have computed the expected ionization rate (which is directly proportional to the doze, which is a quantity responsible for radiation) for different depths of atmosphere (or different altitudes)

In short, the expected ionization rate I is equal to:

![säteilymalli](https://github.com/pyrybjork/Deus-Radiorum/assets/63731201/e3658713-7575-41c4-8520-17b20c287783)

where J(E) is a differential spectrum of cosmic ray particles, Y(E,h) is the so-called ionization yield function, and the summation is over the number of cosmic-ray species.

Since it can look and sound a little bit complicated, I attach a script that allows you to calculate the expected ionization rate (or dose) as a function of atmospheric depth.

For cosmic rays, I took some values from 2021, so it should be a more or less realistic result.

Jos jotakuta kiinnostaa, niin ne käytti siis jotain tämmöstä mallia https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2006JD007150 ja niiden mittausdataa.

## Muonien esiintyvyys

Tosiaan tuo malli ottaa huomioon kaiken säteilyn ja meidän systeemit luultavasti (selvityksen alla) havaitsee vain muoneja. 

Jos näin on malli ei tuli täysin vastaamaan odotettuja mittaustuloksia.

Sergey jakoi kuvaajan, joka näyttää eri hiukkasten suhteellisia määriä eri korkeuksilla.

![CR cascade profile](https://github.com/pyrybjork/Deus-Radiorum/assets/63731201/8f6b2b41-2d10-4280-9c25-09458a1371bd)
