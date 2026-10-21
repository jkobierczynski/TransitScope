// SPDX-License-Identifier: GPL-3.0-or-later
// CFITSIO backend.
//   - image index = n-th HDU that holds a >= 2-D image (empty primaries and tables are skipped,
//     tile-compressed images count), the same numbering libxisfconv uses
//   - BSCALE/BZERO applied by fits_read_img; BLANK pixels become NaN
//   - long strings (CONTINUE convention) are joined into one card
#include <fitsio.h>

#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>

#include "transit/image_reader.h"

namespace transit {

HeaderCard::Type classify_fits_value(const std::string& value);  // reader_factory.cpp

namespace {

void throw_if(int status, const std::string& what) {
    if (!status) return;
    char msg[FLEN_ERRMSG] = {0};
    fits_get_errstatus(status, msg);
    throw std::runtime_error("CFITSIO " + what + ": " + msg);
}

class CfitsioReader final : public ImageReader {
public:
    ~CfitsioReader() override { close(); }

    void open(const std::string& path) override {
        close();
        int status = 0, nhdu = 0;
        fits_open_file(&f_, path.c_str(), READONLY, &status);
        throw_if(status, "open " + path);
        fits_get_num_hdus(f_, &nhdu, &status);
        throw_if(status, "num_hdus");
        for (int h = 1; h <= nhdu; ++h) {
            int type = 0, naxis = 0;
            fits_movabs_hdu(f_, h, &type, &status);
            if (type == IMAGE_HDU) fits_get_img_dim(f_, &naxis, &status);
            throw_if(status, "scan HDUs");
            if (type == IMAGE_HDU && naxis >= 2) images_.push_back(h);
        }
    }
    int hdu_count() const override { return static_cast<int>(images_.size()); }

    std::vector<HeaderCard> header(int img) const override {
        move_to(img);
        int status = 0, nkeys = 0, more = 0;
        fits_get_hdrspace(f_, &nkeys, &more, &status);
        throw_if(status, "hdrspace");
        std::vector<HeaderCard> cards;
        for (int i = 1; i <= nkeys; ++i) {
            char key[FLEN_KEYWORD], val[FLEN_VALUE], com[FLEN_COMMENT];
            fits_read_keyn(f_, i, key, val, com, &status);
            throw_if(status, "read_keyn");
            std::string k = key, v = val;
            if (k == "CONTINUE") continue;                       // folded into the card before it
            const auto amp = v.rfind("&'");
            if (amp != std::string::npos && amp + 2 == v.size()) { // long string: join CONTINUE cards
                char* longstr = nullptr;
                int st2 = 0;
                fits_read_key_longstr(f_, key, &longstr, com, &st2);
                if (!st2 && longstr) {
                    std::string s;
                    for (const char* p = longstr; *p; ++p) { s += *p; if (*p == '\'') s += '\''; }
                    v = "'" + s + "'";
                    fits_free_memory(longstr, &st2);
                }
                int st3 = 0;
                fits_read_keyn(f_, i, key, val, com, &st3);          // restore the position after the lookup
            }
            HeaderCard c{k, v, com, classify_fits_value(v)};
            if (v.empty()) c.type = HeaderCard::Type::Undefined;
            cards.push_back(std::move(c));
        }
        return cards;
    }

    std::vector<long> dims(int img) const override {
        move_to(img);
        int status = 0, naxis = 0;
        long ax[9] = {0};
        fits_get_img_dim(f_, &naxis, &status);
        fits_get_img_size(f_, 9, ax, &status);
        throw_if(status, "img_size");
        std::vector<long> d(ax, ax + naxis);
        while (d.size() > 2 && d.back() == 1) d.pop_back();      // NAXIS3 = 1 is a plain 2-D image
        return d;
    }

    Image pixels(int img) const override {
        auto d = dims(img);
        Image out(static_cast<int>(d[0]), static_cast<int>(d[1]));   // first plane of a cube
        int status = 0, anynul = 0;
        float nul = NAN;
        fits_read_img(f_, TFLOAT, 1, static_cast<LONGLONG>(out.px.size()), &nul, out.px.data(), &anynul, &status);
        throw_if(status, "read_img");
        return out;
    }
    std::string backend() const override { return "cfitsio"; }

private:
    void move_to(int img) const {
        if (img < 0 || img >= static_cast<int>(images_.size())) throw std::out_of_range("CFITSIO: no image " + std::to_string(img));
        int status = 0, type = 0;
        fits_movabs_hdu(f_, images_[img], &type, &status);
        throw_if(status, "movabs_hdu");
    }
    void close() {
        if (f_) { int status = 0; fits_close_file(f_, &status); f_ = nullptr; }
        images_.clear();
    }
    fitsfile* f_ = nullptr;
    std::vector<int> images_;  // 1-based HDU numbers of the image HDUs
};

}  // namespace

std::unique_ptr<ImageReader> make_cfitsio_reader() { return std::make_unique<CfitsioReader>(); }

bool cfitsio_reentrant() { return fits_is_reentrant() != 0; }

void write_fits_float(const std::string& path, const Image& img, const std::vector<HeaderCard>& extra) {
    int status = 0;
    fitsfile* f = nullptr;
    const std::string p = "!" + path;  // overwrite
    fits_create_file(&f, p.c_str(), &status);
    long ax[2] = {img.nx, img.ny};
    fits_create_img(f, FLOAT_IMG, 2, ax, &status);
    fits_write_img(f, TFLOAT, 1, static_cast<LONGLONG>(img.px.size()), const_cast<float*>(img.px.data()), &status);
    for (const auto& c : extra) {
        if (c.type == HeaderCard::Type::String)
            fits_update_key_str(f, c.key.c_str(), c.raw_value.c_str(), c.comment.c_str(), &status);
        else
            fits_update_key_dbl(f, c.key.c_str(), std::stod(c.raw_value), -10, c.comment.c_str(), &status);
    }
    fits_close_file(f, &status);
    throw_if(status, "write " + path);
}

}  // namespace transit
