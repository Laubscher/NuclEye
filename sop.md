# Mise à jour de NuclEye

À lancer dans un terminal, depuis le dossier du dépôt NuclEye.

## 1. Préparer la mise à jour

Il faut Python 3.9+, MAFFT 7 installé, le taxid NCBI du virus et une adresse e-mail.

Remplacez les valeurs entre chevrons. Exemple : `--since 2026-10` récupère les séquences rendues publiques depuis octobre 2026.

```bash
VIRUS=<nom_du_virus>
TAXID=<taxid_NCBI>
EMAIL=<email>
DEPUIS=2026-10
```

## 2. Vérifier puis mettre à jour les séquences

Vérification locale, sans connexion à Genbank :

```bash
python3 scripts/update_genbank.py "data/$VIRUS" --taxid "$TAXID" --email "$EMAIL" --since "$DEPUIS" --check
```

Lancer la mise à jour GenBank :

```bash
python3 scripts/update_genbank.py "data/$VIRUS" --taxid "$TAXID" --email "$EMAIL" --since "$DEPUIS"
```

Optionnel : ajouter `--dry-run` à la commande pour rechercher et contrôler les séquences sans modifier les fichiers.

## 3. Régénérer les JSON

```bash
python3 scripts/update_data.py
```

## 4. Publier la mise à jour


Ajouter les JSON utilisés par la page, puis pousser sur `main` :

```bash
git add data/viruses.json "data/$VIRUS"/*/data.json
git commit -m "Update NuclEye data"
git push origin main
```

GitHub Actions republie ensuite la page.
