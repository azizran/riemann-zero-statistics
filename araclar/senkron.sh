#!/bin/bash
# senkron.sh — iki Mac (MacBook + Mac mini) arası senkron.
#   ./araclar/senkron.sh           oturum başı: GitHub'dan çek + köprüleri kur + durum
#   ./araclar/senkron.sh --gonder  oturum sonu: yukarıdakiler + yerel commit'leri push
# Kod/defter GitHub üzerinden, veri (ara önbellekler) YENİDEN ÜRETİMLE senkronlanır
# (araclar/yeniden_kur.sh). Hiçbir git-takipli dosyayı değiştirmez; force-push YOK.
set -u

REPO="$(cd "$(dirname "$0")/.." && pwd)"
KALICI="$REPO/qm_riemann/scratchpad"
UUID="71b922d7-9b7c-485d-9927-3b0db425e2a2"
ESKI_TMP="/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/$UUID/scratchpad"
ESKI_HOME="/Users/ugursezen/Desktop/arin/deney"
GONDER=0; [ "${1:-}" = "--gonder" ] && GONDER=1

echo "== makine: $(scutil --get ComputerName 2>/dev/null || hostname)   repo: $REPO"

echo; echo "== 1) git =="
cd "$REPO" || exit 1
git fetch -q origin || { echo "  UYARI: fetch başarısız (internet?)"; }
KIRLI=$(git status --porcelain --untracked-files=no | wc -l | tr -d ' ')
read ILERI GERI < <(git rev-list --left-right --count HEAD...@{u} 2>/dev/null || echo "0 0")
echo "  yerel ileride: $ILERI   geride: $GERI   commit'lenmemiş değişiklik: $KIRLI"
if [ "$ILERI" -gt 0 ] && [ "$GERI" -gt 0 ]; then
  echo "  !! AYRIŞMA: iki makine farklı commit'ler üretmiş. Otomatik birleştirme YAPILMADI."
  echo "     Çözüm (Claude'a sor ya da): git merge origin/main  → sonra --gonder"
elif [ "$GERI" -gt 0 ]; then
  if [ "$KIRLI" -gt 0 ]; then
    echo "  !! Geride ama commit'lenmemiş değişiklik var — çekme YAPILMADI (önce commit/stash)."
  else
    git merge -q --ff-only origin/main && echo "  çekildi: $GERI commit (fast-forward)."
  fi
fi
if [ "$ILERI" -gt 0 ] && [ "$GERI" -eq 0 ]; then
  if [ "$GONDER" -eq 1 ]; then git push -q origin main && echo "  gönderildi: $ILERI commit."
  else echo "  not: $ILERI commit gönderilmedi (oturum sonunda: senkron.sh --gonder)."; fi
fi
[ "$ILERI" -eq 0 ] && [ "$GERI" -eq 0 ] && echo "  GitHub ile senkron."

echo; echo "== 2) veri köprüsü (eski geçici yol → kalıcı klasör) =="
mkdir -p "$KALICI"
if [ -L "$ESKI_TMP" ]; then
  ln -sfn "$KALICI" "$ESKI_TMP"; echo "  köprü var: $ESKI_TMP -> $KALICI"
elif [ -d "$ESKI_TMP" ]; then
  # gerçek klasör: önce kalıcıya kopyala, farksızsa köprüye çevir
  rsync -a "$ESKI_TMP/" "$KALICI/"
  FARK=$(rsync -ain "$ESKI_TMP/" "$KALICI/" | wc -l | tr -d ' ')
  if [ "$FARK" -eq 0 ]; then
    rm -rf "$ESKI_TMP" && ln -s "$KALICI" "$ESKI_TMP"
    echo "  geçici klasör kalıcıya taşındı ve köprüye çevrildi."
  else
    echo "  !! kopyada $FARK fark var — geçici klasöre DOKUNULMADI."
  fi
else
  mkdir -p "$(dirname "$ESKI_TMP")" && ln -s "$KALICI" "$ESKI_TMP"
  echo "  köprü kuruldu (sistem geçici klasörü silmişti): $ESKI_TMP -> $KALICI"
fi

echo; echo "== 3) eski betiklerin mutlak yolu =="
if [ -e "$ESKI_HOME/qm_riemann" ]; then
  echo "  $ESKI_HOME erişilebilir — TAMAM."
else
  echo "  Bu makinede eski yol yok (Mac mini). BİR KEZ şunu çalıştır (yönetici şifresi ister):"
  echo "    sudo mkdir -p /Users/ugursezen/Desktop/arin"
  echo "    sudo ln -sfn \"$REPO\" $ESKI_HOME"
fi

echo; echo "== 4) veri durumu =="
EKSIK=0
for f in 155/eta_son_t0.4_c4000.npz 155/z_Hkeskin.npy 184/K1_gercek.npz 185/K1_faktorler.npz \
         186/G1_proj_gercek.npz 187/K2_zeta.json 188/harita_K_gercek.npz; do
  if [ -e "$KALICI/$f" ]; then echo "  var   $f"; else echo "  YOK   $f"; EKSIK=1; fi
done
[ "$EKSIK" -eq 1 ] && echo "  → eksikler için: ./araclar/yeniden_kur.sh  (~40 dk; 187 sayılarını kapıyla doğrular)"
du -sh "$KALICI" 2>/dev/null | awk '{print "  kalıcı veri: " $1}'
